#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}


const double degreeFactor(-deg2rad(25));

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;



  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
//    Disable input cout
//    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double throttle = j[1]["throttle"];
//          double psi_u = j[1]["psi_unity"];
          double steer = j[1]["steering_angle"];
          double v = j[1]["speed"];


          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value;
          double throttle_value;

          //first step: convert all coordinates into car coords
          //Therefore we must use the invers of the rotation matrix R-1
          //This can either be done by changing the Rotation matrix or
          //simply by using the negative angle for rotation R-1(alpha) = R(-alpha)
          //refer to https://de.wikipedia.org/wiki/Drehmatrix e.g.
          //Furthermore the car coord system starts at 0,0. That means
          //we'll only identify the difference from px,py to the waypoints
          //==> dx, dy and start from 0,0 as point
          double const car_x(0.), car_y(0.);
          vector<double> car_ptsx(ptsx.size(), 0.);
          vector<double> car_ptsy(ptsy.size(), 0.);
          //the car has no heading - spi 0i01s zero
          //TODO: previously it was zero
//          double const car_psi(steer*degreeFactor);
          double const car_psi(0.);

          assert(ptsx.size()==ptsy.size());
          //use negative psi to the the rotation matrix invers
          double const sinPsi(sin(-psi));
          double const cosPsi(cos(-psi));
          //for further processing, we immediately convert the resulting
          //points into an eigen-representation
//          Eigen::VectorXd ePtsX(car_ptsx.size()), ePtsY(car_ptsy.size());
          vector<double> car_X, car_Y;
          for(size_t i(0U), _maxI(ptsx.size()); i<_maxI;i++)
          {
            //what is the delta from waypoints to car-coord
            double const dpx(ptsx[i]-px);
            double const dpy(ptsy[i]-py);
            //apply rotation matrix
            car_ptsx[i] = ( car_x + (cosPsi*dpx - sinPsi*dpy));
            car_ptsy[i] = ( car_y + (sinPsi*dpx + cosPsi*dpy));
            if(/* (car_ptsx[i]>=0.) &&*/ (car_ptsx[i]<80.))
            {
              car_X.push_back(car_ptsx[i]);
              car_Y.push_back(car_ptsy[i]);
            }
            //eigen represetnation - only consider points in front of the car
//            if((car_ptsx[i] >= 0.))// && (car_ptsy[i] >= 0.))
//            {
//              ePtsX[i] = car_ptsx[i];
//              ePtsY[i] = car_ptsy[i];
//            }
          }
          Eigen::VectorXd ePtsX(car_X.size()), ePtsY(car_X.size());
          for(int i(0), _maxI(car_X.size()); i < _maxI; i++)
          {
            ePtsX[i]=car_X[i];
            ePtsY[i]=car_Y[i];
          }
          cout<<"SIZE "<<ePtsX.size()<<endl;
          //we are searching for a 2nd order polynomial
          auto coeffs = polyfit(ePtsX, ePtsY, 2);
          assert(coeffs.size()==3);
          cout <<" Coeffs "<<coeffs<<endl;

          //calculate the CTE
          //  = distance of f(x) and y
//          for(int i(0); i < car_ptsx.size(); i++)
//          {
//            testo.push_back(polyeval(coeffs, car_ptsx[i]));
//          }

          //calculate the orientation error according to chapter 9
          //take care of 2nd order polynom
          double const epsi(car_psi - atan(coeffs[1]/*+2*coeffs[2]*px -> we want to know the angle at the car-position - so x is more or less zero, we are only interested in the ypart*/));
          double const cte(polyeval(coeffs, car_x) - car_y);

          double const angSpeed(v/2.67); //angular velocity

          double const delta_psi(v / 2.67 * -steer * 0.1); //see 18-4
          double const delta_epsi(delta_psi - atan(coeffs[1]));
          double const delta_v(throttle * 0.1);
          double const delta_cte(sin(epsi) * v * 0.1);
          double const delta_x = v * cos(-steer) * 0.1;
          double const delta_y = v* sin(-steer) * 0.1;
          cout<<"CTE "<<cte <<" DeltaCTE "<<delta_cte<<" DX "<<ePtsX[0] << " DY "<<ePtsY[0]<<endl;


          //fill stateVector
          Eigen::VectorXd state(6);
          state << car_x+delta_x, car_y+delta_y, car_psi+delta_psi, v+delta_v, cte+delta_cte, epsi+delta_epsi;

//          cout << "State >> x ,y ,psi ,v : "<<car_x<<", "<<car_y<<", "<<car_psi<<", "<<v<<endl;
//          cout << "\t cte, epsi :"<<cte<<", "<<epsi<<endl;

          //do solving
          auto vars = mpc.Solve(state, coeffs);

          //Instead of changing the the equation for the steering input
          //refer to Project- Tricks, we simply change the steering value
          double const nominator(deg2rad(25));
          cout<<"SteeringValue: "<<vars[0]<<endl<<endl;
          steer_value=-1*vars[0]/nominator;
//          if(v<10.)
//          {
            throttle_value=vars[1];
//          }
//          throttle_value=vars[1];
//          cout<< "SteeringV "<<vars[0]<<" Throttle "<<vars[1]<<endl;


          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          for(int i(2); i<vars.size(); i+=2)
          {
            mpc_x_vals.push_back(vars[i]);
            mpc_y_vals.push_back(vars[i+1]);
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          next_x_vals = car_X;//car_ptsx;
          next_y_vals = car_Y;//car_ptsy;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          //LATENCY consider this post https://discussions.udacity.com/t/calibration-for-the-acceleration-and-steering-angle-for-latency-consideration/276413

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
//          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
