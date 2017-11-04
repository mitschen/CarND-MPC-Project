# CarND-MPC Control Project


###Overview
This project is about the implementation of a MPC (model predicted controller) which interacts with the Udacity car simulator. The MPC should provide steering angle and acceleration based on the state information provided by the simulator. Goal of the MPC is to make sure that the car stays on track.

### Implementation
The implementation is splitted into two files. The ``main.cpp`` contains the interface (websocket) to the car simulator. Is gathers the state information, transforms them into a car- state and forwards it to the MPC. The result of the MPC is converted into the car-simulator units and send back to the simulator. In order to simulate process delay, a artificial sleep of 0.1 second is used before the result is send back to the car simulator.

The other file is the ``mpc.cpp`` which represents the MPC itself. Based on the input state it calculates a throttle- and steering value output, trying to optimize the error-cost given by different dynamic-equations.

The next chapter will go through the different implementation parts in more details.

#### The Model
The start of the dataflow is the input, provided by the car-simulator: 
- the waypoint list as global positions (map origin) as x,y coordinates: **ptsx, ptsy**
- the orientation of the car with respect to the map origin: **psi**
- the global coordinates of the vehicle: **px, py**
- the current steering angle of the car in radians: **steer**
- the throttle of the car: **throttle**
- and the current speed of the car: **v**

Due to the fact that all our calculations and projections should happen in car-coordinate system, we start with the first step converting the state from map into car coordinate system.
In car coordinate system we're assuming that x and y coordinates are zero and the psi (orientation of the car) is zero as well.

```
double const car_x(0.), car_y(0.);
double car_psi(0.);
```

But - in order to project the waypoint list in the car coordsystem we need to transform them fron map origin to car origin. Therefore we calculate the difference of the waypoints with respect to the car position and apply an inverse rotation of the coordinates to the car origin.

```
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
```

With this transformation we get a waypoint coordinate list based on the car-coordinates origin. The ``car_X, car_Y`` values are used later to draw the yellow reference line in front of the car.

Given the waypoints in car-coords, we use the ``polyfit`` helper function to identify a suitable polygon for the waypoint trajectory. My first setup was using a 2nd order polynom but i was facing some weird results from the solver. Furthermore the solver trajectory was not matching very well with the waypoints. Order 3 polynom did fit much better. The result of the ``polyfit`` function is a list of 4 coefficients. These coefficients are needed to predict the dynamic behaviour of the car in the MPC later on.

#####latency consideration
Instead of feeding the MPC with the status quo given by the simulator, I already did the prediction of the car position considering the process delay of 100ms. So I'm calculating the car-position ``car_x, car_y`` as well as the new car orientation ``car_psi`` based on steering angle, throttle, velocity and delay. 
```
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
```

With this predicted future state the ``cte`` and the ``error-psi`` is given according to the following equation:
```
double const epsi(car_psi - atan(deriviate(coeffs, delta_x)));
double const cte(polyeval(coeffs, delta_x) - delta_y);
```

The final state vector that is used as input for the MPC is therefore
``state << delta_x, delta_y, car_psi, car_v, cte, epsi``


####MPC processing
The MPC is more or less following the suggestion given in Chapter 19. The MPC tries to optimize the error-value under consideration of a certain cost-function and a certain time span.
For the timespan I choose ``c_N=20`` steps with a cycle time of ``c_dt=0.1`` seconds - see ``MPC.h``. That means in other words the MPC tries to minimize the overall cost value for a timeperiod of 2 seconds in future. Please note: the process latency is already normalized in previous step.

The solver is considering the dynamic car model that was discussed in chapter 12 (simple bike model). Based on the states at t it optimizes states at t+1 considering the dynamic model and the cost function.

```
if(yawrate is zero):
	x_[t+1] = x[t] + v[t] * cos(psi[t]) * c_dt
	y_[t+1] = y[t] + v[t] * sin(psi[t]) * c_dt
if(yawrate != zero):
	x_[t+1] = x_[t] + v[t] / c_Lf * delta[t] * (sin(psi[t] + v[t] / c_Lf * delta[t] * c_dt) - sin(psi[t]))
	y_[t+1] = y_[t] + v[t] / c_Lf * delta[t] * (-cos(psi[t] + v[t] / c_Lf * delta[t] * c_dt) + cos(psi[t]))
psi_[t+1] = psi[t] + v[t] / c_Lf * delta[t] * c_dt
v_[t+1] = v[t] + a[t] * c_dt
cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * c_dt
epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / c_Lf * c_dt
```

As cost function I stay with the suggestion from chapter 19 adding some cost-optimizations:
The costvalue consists of differen parts:
In the first part, the absolute difference of y, psi and speed (=cte, epsi, speed-diff) is taken into consideration. I figured out very fast, that the speed-impact must be reduced in order to get a suitable result from the solver. In other words, optimizing the speed is not that important than optimizing the cte or epsi error.
Furthermore i've adjusted the weights for the cte and epsi. My consideration here was, that my car should not hit the optimal matching with the waypoint trajectory at time t, but it should be optimal in the very close future. In other words instead of using the same weight for all timesteps, I increased the penalty-costs for the close future, reduced it for current state and states in far future. In sourcecode this looks like:
```
fg[0] += (CppAD::pow(vars[cte_start + t], 2)) / (1.+ pow(double(t-5), 2));
fg[0] += (CppAD::pow(vars[epsi_start + t], 2)) / (1.+ pow(double(t-5), 2));
//reduction of speed-weight - otherwise the solver provides really
//weird curvatures
fg[0] += 0.1 * CppAD::pow(c_ref_v - vars[v_start + t] , 2) ;
```
The nominator ``(1. +pow(double(t-5), 2))`` makes sure, that the timestep t+5 has the biggest impact on the cost-value. That means futhermore, the solver tries to optimize this timestep which results in the position 0.5 seconds in future more then any other one.

The second part is the unchanged minimization of the actuators - meaning resulting steering and throttle.
```
fg[0] += CppAD::pow(vars[delta_start + t], 2);
fg[0] += CppAD::pow(vars[a_start + t], 2);
``` 
And finally the third part considers the change of throttle and steering angle between two timestamps. 
```
fg[0] += 200.0*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
```

I was facing a heavy thrashing of the steering angle with the original weights, so i've increase the weight of the steering angle diffs by 200 to achieve a smoother behaviour.


###loopback to simulator
As a result of the MPC we are getting the steering value as well as the throttle values, combined with the predicted x,y coordinates of the car. The steering value is normalized by the constant 25 degree to gurantee a steering value between [-1, 1]. Throttle is forwarded unchanged and the trajectory is drawn as a green line in the simulator. 


### Summary
The project realizes a MPC which steers a car through a track. The MPC uses motion model and current state parameters to predict the car movement in future. Based on this prediction proper steering and throttle values are derived which will be send back to the simulation engine. 
The artificial processing delay of 100ms was compensated by the prediction of the car position. This predicted car position was used as initial state parameter for the MPC.

The MPC is a controller with hundrets of hyperparameters. For me it wasn't easy to tweak this MPC so that it runs the way I want it to go. The major problem I'm facing is the solver and the resulting trajectory. Sometimes the trajectory is resulting in circles, sometimes it is pointing in the complete opposite direction of the waypoints-trajectory. Unfortuantely I still do not understand this behaviour of the solver. Combined with the fact, that the simulator doesn't allow any debugging attempts - I've overcome the problem only by try and error. Still I didn't explore the impact of each hyperparameter, like reducing  the number of steps ``c_N`` or changing the cost function.
