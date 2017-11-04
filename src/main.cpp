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
double polyeval(Eigen::VectorXd coeffs, double x)
{
  double result(0.0);
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}
// Evaluate a derivate of the polynomial.
double deriviate(Eigen::VectorXd coeffs, double x)
{
  double result(0.);
  for(int i(1); i<coeffs.size(); i++)
  {
    result += coeffs[i] * pow(x, (i-1)) * i;
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
          double psi_u = j[1]["psi_unity"];
          double steer = j[1]["steering_angle"];
          double v = j[1]["speed"];
          //no used variable
          (void) psi_u;

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value(0.);
          double throttle_value(0.);

          //first step: convert all coordinates into car coords
          //Therefore we must use the invers of the rotation matrix R-1
          //This can either be done by changing the Rotation matrix or
          //simply by using the negative angle for rotation R-1(alpha) = R(-alpha)
          //refer to https://de.wikipedia.org/wiki/Drehmatrix e.g.
          //Furthermore the car coord system starts at 0,0. That means
          //we'll only identify the difference from px,py to the waypoints
          //==> dx, dy and start from 0,0 as point
          double const car_x(0.), car_y(0.);

          //the car has no heading - spi 0i01s zero
          double car_psi(0.);

          assert(ptsx.size()==ptsy.size());
          //use negative psi to the the rotation matrix invers
          double const sinPsi(sin(-psi));
          double const cosPsi(cos(-psi));
          //for further processing, we immediately convert the resulting
          //points into an eigen-representation
          vector<double> car_X(ptsx.size(), 0.), car_Y(ptsy.size(), 0.);
          for(size_t i(0U), _maxI(ptsx.size()); i<_maxI;i++)
          {
            //what is the delta from waypoints to car-coord
            double const dpx(ptsx[i]-px);
            double const dpy(ptsy[i]-py);
            //apply rotation matrix
            car_X[i] = ( car_x + (cosPsi*dpx - sinPsi*dpy));
            car_Y[i] = ( car_y + (sinPsi*dpx + cosPsi*dpy));
          }
          //translate into eigen representation
          Eigen::VectorXd ePtsX(car_X.size()), ePtsY(car_X.size());
          for(int i(0), _maxI(car_X.size()); i < _maxI; i++)
          {
            ePtsX[i]=car_X[i];
            ePtsY[i]=car_Y[i];
          }
          const int order_poly(3);

          //second step: calculate the polynom which represents the
          //reference path in car-coordinates
          auto coeffs = polyfit(ePtsX, ePtsY, order_poly);
          assert(coeffs.size()==order_poly+1);


          //third step: forecast the position, error,... for future
          //assuming a delay of 0.1 seconds
          //see/compare the equation from lesson 12-3 and 18-4
          double const yawrate(v * -steer / c_Lf);
          double car_v(v+throttle * c_delay);//see 18-4

          //Delta_xyz is the offset on the car_xyz position
          double delta_x = 0.;
          double delta_y = 0.;
          //if yaw rate is zero we have to consider another equation
          //Due to the fact that this method calculates the delta under
          //consideration of time and yawrate, we're using the car_psi which
          //is zero
          if(pow(yawrate, 2.)>0.00001)
          {
            delta_x = v / yawrate * ( sin(car_psi + yawrate * c_delay) - sin(car_psi) );
            delta_y = v / yawrate * ( cos(car_psi) - cos(car_psi + yawrate * c_delay) );
          }
          else
          {
            delta_x = v * 1.0 * c_delay; // v * cos(0) * delay
            delta_y = 0.; // v * sin(0) * delay
          }
          //now calculate the car_psi based on yawrate and the delay
          car_psi += yawrate * c_delay;

          //calculate the orientation error according to chapter 9
          //Please note: we are using delat_x/delta_y due to the fact, that
          //car_x, car_y is zero and we want to predict the future error
          double const epsi(car_psi - atan(deriviate(coeffs, delta_x)));
          double const cte(polyeval(coeffs, delta_x) - delta_y);

          //fill stateVector
          Eigen::VectorXd state(6);
          state << car_x+delta_x, car_y+delta_y, car_psi, car_v, cte, epsi;

#if 0 //debuggin purpose only
          cout << "State >> x ,y ,psi ,v : "<<car_x+delta_x<<", "<<car_y+delta_y<<", "<<car_psi<<", "<<car_v<<endl;
          cout << "\t cte, epsi :"<<cte<<", "<<epsi<<endl;
#endif
          //do solving
          auto vars = mpc.Solve(state, coeffs);

          double const nominator(deg2rad(25));
          //Instead of changing the the equation for the steering input
          //refer to Project- Tricks, we simply change the steering value
          steer_value=-1*vars[0]/nominator;
          throttle_value=vars[1];


          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          for(auto i(2u); i<vars.size(); i+=2u)
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
          this_thread::sleep_for(chrono::milliseconds((int)(c_delay * 1000)));
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
