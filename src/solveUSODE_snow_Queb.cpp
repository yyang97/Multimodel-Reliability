#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <stdlib.h>
#include <math.h>
#include <iostream> // for cin and cout
#include <fstream>
#include <sstream>  
#include <vector> // for vector
#include <map>
#include <string>
#include <omp.h>
#include <random> // for distributions
#include "boost/math/special_functions/gamma.hpp" // for the incomplete gamma function
#include "boost/numeric/odeint.hpp" // for the ode solver
#include "boost/math/tools/roots.hpp" // for the root finding

using namespace std;
using namespace boost::numeric::odeint;
using boost::math::tgamma_lower;
using namespace boost::math::tools; // for bisect

////////////////////////////////////////////////////////////////////////////////

// setup parameters
const double tau0 = 6467.75; // median of sample tau_s
const double mu = 1;
const double k = 388440;
const double T0 = 4500 / k;
const int n_per_theta = 100000;
const eps_tolerance<double> tol(16);
const double t_start = 0;
const double t_end = 438300; // hours in 50 years
double dt = 100; // integration step size

map<string, double> load_profile_par{
  {"phi", 1},
  {"load_s_shape", 3.122},
  {"load_s_scale", 0.0481},
  {"load_p_shape", 0.826},
  {"load_p_scale", 0.1023},
  {"R0", 3000},
  {"alpha_d", 1.25},
  {"alpha_l", 1.5},
  {"load_d_mean", 1.05},
  {"load_d_sd", 0.1},
  {"gamma", 0.25},
  {"mean_Ts", 10},
  {"mean_Te", 0.5833},
  {"mean_Tp", 0.03835},
  {"N", 100},
  // Here are the parameters for snow
  {"NS", 10},
  
  // ground snow load parameters
  // Quebec
  {"A", 0.3222}, 
  {"B", 17.0689},		
  // ground-to-roof transformation
  {"r_flat_mean", 0.6},
  {"r_flat_COV", 0.45}  
  
  };

// set up random number generators
mt19937 generator(0);

////////////////////////////////////////////////////////////////////////////////

double stepfun(double x, vector<double> x_vec, vector<double> y_vec) {
  if (y_vec.size() - x_vec.size() != 1) {
    cout << "The length of y_vec must be greater than the length of x_vec by 1." << endl;
    return (0);
  }
  int n = x_vec.size();
  if (x < x_vec[0]) return y_vec[0];
  else if (x >= x_vec[n - 1]) return y_vec[n];
  else {
    for (int i = 1; i < n; ++i) {
      if (x < x_vec[i] && x >= x_vec[i - 1]) return y_vec[i];
    }
  }
}

double sigmalog (double mu, double sigma){
  double result = sqrt(log(1 + exp(log(sigma*sigma) - 2 * log(mu))));
  return result;
}


class load_profile {
  double load_p_shape, load_p_scale, load_s_shape, load_s_scale, phi, R0, alpha_d, alpha_l,
  mean_Ts, mean_Te, mean_Tp, D_d, load_d_mean, load_d_sd, gamma, A, B,
  r_flat_mean,r_flat_COV;
  int N = 550;
  int NS = 10;
  
  // probability of snow in a segment.
  double p_e =  1.00 - exp(-1.00/NS*(exp(A*B)));  
  double p_0 = exp(-exp(A*B)); // probability of no snow in one year
  
  //double Bstr = A*B/(A*B+3.3843);
  //double Astr = A*B+3.3843;
  
  double Bstr = A*B/(A*B+3.9019); // parameters for distribution g
  double Astr = A*B+3.9019;

public:
  vector<double> T_s = vector<double>(N);
  vector<double> T_e = vector<double>(N);
  vector<double> load_s = vector<double>(N + 1);
  vector<double> load_e = vector<double>(N + 1);

  load_profile(map<string, double> load_profile_par) :
  load_s_shape(load_profile_par["load_s_shape"]),
  load_s_scale(load_profile_par["load_s_scale"]),
  load_p_shape(load_profile_par["load_p_shape"]),
  load_p_scale(load_profile_par["load_p_scale"]),
  phi(load_profile_par["phi"]),
  alpha_d(load_profile_par["alpha_d"]),
  alpha_l(load_profile_par["alpha_l"]),
  R0(load_profile_par["R0"]),
  load_d_mean(load_profile_par["load_d_mean"]),
  load_d_sd(load_profile_par["load_d_sd"]),
  gamma(load_profile_par["gamma"]),
  mean_Ts(load_profile_par["mean_Ts"]),
  mean_Te(load_profile_par["mean_Te"]),
  mean_Tp(load_profile_par["mean_Tp"]),
  A(load_profile_par["A"]),
  B(load_profile_par["B"]),
  r_flat_mean(load_profile_par["r_flat_mean"]),
  r_flat_COV(load_profile_par["r_flat_COV"])
    
  
  
  
  {
    normal_distribution<double> norm_d(load_d_mean, load_d_sd);
    gamma_distribution<double> gamma_s(load_s_shape, load_s_scale);
    gamma_distribution<double> gamma_p(load_p_shape, load_p_scale);
    exponential_distribution<double> exp_Ts(1.0 / mean_Ts);
    exponential_distribution<double> exp_Te(1.0 / mean_Te);
    exponential_distribution<double> exp_Tp(1.0 / mean_Tp);
    
    // snow simulation
    uniform_real_distribution<double> runif(0,1);
    double mu = r_flat_mean;
    double sigma = r_flat_COV * r_flat_mean;
    double input1 = log(mu) - sigmalog(mu, sigma)*sigmalog(mu, sigma) / 2;
    double input2 = sigmalog(mu, sigma);
    lognormal_distribution<double>log_norm_r(input1,input2);
    
    
    
    
    
    
    // generate the dead load
    D_d = norm_d(generator);
    // generate the non-winter time
    T_e[0] = exp_Te(generator);
    // there is no snow, snow load = 0
    load_e[0] = 0;
    
    for(int i = 1; i < N; ++i){
      if(i % 11 != 0){
        // there  might be snow in this segment as it is winter
        double r_n = runif(generator);
        // generate the time
        T_e[i] = T_e[i-1] + exp_Tp(generator);				  
        if(r_n > p_e){
          // there is no snow in that segment
          load_e[i] = 0;
        }
        else{
          // there is snow in this segment
          double p = runif(generator);
          double r = log_norm_r(generator);
          load_e[i] = (Bstr + 1.0 / Astr * (-log(-NS * log(1.0 - p_e + p_e * p))))*r;			    
        }
      }
      else{
        // there is no snow
        T_e[i] = T_e[i-1] + exp_Te(generator);
        load_e[i] = 0;
      }
    }
    
    // Snow in the last segment is 0.
    load_e[N] = 0;
  
  
  
  
  }

  double operator()(double t) {
    t = t / 8760; // express t in years
    double D_e = stepfun(t, T_e, load_e);
    double load_t = phi * R0 * (gamma * D_d + D_e)/(gamma * alpha_d + alpha_l);
    return load_t;
  }

  vector<double> max() {
    double max_load;
    vector<double> res(2, 0);
    for (double t = t_start; t < t_end; t += dt) {
      max_load = (*this)(t);
      if (max_load > res[1]) {
        res[0] = t;
        res[1] = max_load;
      }
    }
    return (res);
  }

  void write(ofstream& file){
    if(file.is_open()){
      file << "T_e, load_e" << endl;
      
      for(int i = 0; i < T_s.size(); ++i){
        file << T_e[i] << "," << load_e[i] << endl;
        
      }	
      file << " ," << load_e[T_s.size()] << endl;
      
    }
    else{
      cout << "File is not open!" << endl;
    }
  }

};

class constant_load {
  double load_level;

public:

  constant_load(double l) : load_level(l) {
  };

  double operator()(double t) {
    // t is in hour
    double T0 = load_level / k;
    if (t <= T0) return k * t;
    else return load_level;
  }

  vector<double> max() {
    double T0 = load_level / k;
    vector<double> res(2, 0);
    res[0] = T0;
    res[1] = load_level;
    return (res);
  }

};

template <class T>
void USADM(const double &x, double &dxdt, const double t, vector<double> parameter, T tau) {

  double a, b, tau_s;

  a = parameter[0];
  b = parameter[1];   // B' = B / tau0
  tau_s = parameter[2];  // also already divided by tau0

  dxdt = exp( -a + b * tau(t) / tau_s);
  
}

class BadConversion : public std::runtime_error {
public:

  BadConversion(const std::string& s)
  : std::runtime_error(s) {
  }
};

double convertToDouble(const string& s) {
  istringstream i(s);
  double x;
  char c;
  if (!(i >> x))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return x;
}

vector< vector<double> > readCSV_to_vec(ifstream& file) {
  vector< vector<double> > result;
  vector<double> row;
  double p;
  string line;
  string cell;

  while (getline(file, line)) {

    stringstream lineStream(line);

    while (getline(lineStream, cell, ',')) {
      p = convertToDouble(cell);
      row.push_back(p);
    }
    result.push_back(row);
    row.clear();
  }
  return result;
}

vector<double> generate(vector<double> &theta) {
  // Generate random variates
  normal_distribution<double> normZ(0, 1);

  vector<double> pars(3);
  pars[0] = theta[0]; // A
  pars[1] = theta[1]; // B' = B / tau_0
  pars[2] = exp( theta[2]* normZ(generator));

  return pars;
}



////////////////////////////////////////////////////////////////////////////////


// main function

int main(int argc, char* argv[]) {

  double x, t;
  double x_max = 1;


  if (argc == 3) load_profile_par["phi"] = atof(argv[2]);
  cout << "phi is set to be " << load_profile_par["phi"] << endl;

  //cout << tol << endl;  

  // read theta from file
  ifstream fin(argv[1]);
  if (fin.is_open()) {
    vector< vector<double> > theta = readCSV_to_vec(fin);
    int ntheta = theta.size();

    vector<int> nFail_per_theta(ntheta, 0);
    vector<int> nFail_per_theta_noDOL(ntheta, 0);
    vector<double> prop_Fail(ntheta, 0);
    vector<double> prop_Fail_noDOL(ntheta, 0);
    vector<double> time_to_failure(ntheta * n_per_theta);

    vector<double> pars;

    cout << "Read " << ntheta << " thetas!" << endl;
    cout << "Start simulation!" << endl;



    omp_set_num_threads(16);
#pragma omp parallel for schedule(dynamic) private(t, x, pars)  
    for (int i = 0; i < ntheta; ++i) {
      int tmp = 0;
      int tmp_noDOL = 0;


      adams_bashforth_moulton< 5, double > stepper;

      for (int j = 0; j < n_per_theta; ++j) {



        pars = generate(theta[i]);

        load_profile tau(load_profile_par);

        vector<double> max_t = tau.max();

          // no DOL
          if (max_t[1] > pars[2] * tau0) {
            tmp_noDOL += 1;
            tmp += 1;
            continue;
          }

          // solve the ODE 		
          auto model = [&pars, &tau](const double &x, double &dxdt, const double t) {
            USADM(x, dxdt, t, pars, tau);
          };


          x = 0; // initial state
          t = t_start;
          //times.push_back(t);
          //x_vec.push_back(x);
          stepper.initialize(model, x, t, dt);
          t += dt;
          while (t < t_end) {
            stepper.do_step(model, x, t, dt);
            //times.push_back(t);
            //x_vec.push_back(x);
            if (x > x_max || std::isnan(x)) // if x > 1, stop calculation
            {
              tmp += 1;
              break;
            }
            t += dt;
          }

          time_to_failure[i * n_per_theta + j] = t / 8760.0;


        

        // the ode can be solved simply with the following; but we want it stops early when some criterion is met
        // so we have to solve it manually.

        //integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) ,
        //    model , x , t_start , t_end , dt, push_back_state_and_time(x_vec, times));

      }
      nFail_per_theta[i] = tmp;
      nFail_per_theta_noDOL[i] = tmp_noDOL;
      prop_Fail[i] = (double) nFail_per_theta[i] / (double) n_per_theta;
      prop_Fail_noDOL[i] = (double) nFail_per_theta_noDOL[i] / (double) n_per_theta;

      cout << i << " " << tmp << " " << tmp_noDOL << endl;

    }

    // the output file names depend on phi
    stringstream prob_name;
    stringstream prob_noDOL_name;
    stringstream time_name;

    prob_name << "USprobQueb_" << load_profile_par["phi"] << ".csv";
    prob_noDOL_name << "USprobQueb_noDOL_" << load_profile_par["phi"] << ".csv";
    time_name << "UStimeQueb_" << load_profile_par["phi"] << ".csv";


    // output time_to_failure
    ofstream file2(time_name.str());
    if (file2.is_open()) {
      for (int i = 0; i < time_to_failure.size(); ++i) {
        file2 << time_to_failure[i] << endl;
      }
    } else {
      cout << "File is not open!" << endl;
    }


    ofstream file3(prob_name.str());
    if (file3.is_open()) {
      for (int i = 0; i < prop_Fail.size(); ++i) {
        file3 << prop_Fail[i] << endl;
      }
    } else {
      cout << "File is not open!" << endl;
    }


    ofstream file4(prob_noDOL_name.str());
    if (file4.is_open()) {
      for (int i = 0; i < prop_Fail.size(); ++i) {
        file4 << prop_Fail_noDOL[i] << endl;
      }
    } else {
      cout << "File is not open!" << endl;
    }


  } else {
    cout << "Please specify the file of theta!" << endl;
  }

  return (0);
}

////////////////////////////////////////////////////////////////////////////////



