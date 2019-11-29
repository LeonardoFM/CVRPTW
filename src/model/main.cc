#include <string>
#include <tuple>
#include <vector>
#include <algorithm>
#include <numeric>
#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_enums.pb.h"
#include "ortools/constraint_solver/routing_index_manager.h"
#include "ortools/constraint_solver/routing_parameters.h"

namespace operations_research {
  class DataModel {
  public:
    // the set store all model data in
    const int time_horizon = 24*60*60;
    const std::vector<std::vector<int64>> time_matrix{
      // 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 16
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //0
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //1
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //2
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //3
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //4
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //5
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //6
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //7
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //8
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //9
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //10
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //11
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, //12
        {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}}; //13        
                                 //    0  1  2  3  4  5  6  7  8   9   10  11 12 13 14 15 16
      const std::vector<int64> demands{0, 0, 0,15,15,15, 25,25,25,-15,-15,-15, 0, 0};
      const std::vector<int64> demand_load_time{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      const std::vector<int64> demand_unload_time{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      const std::vector<int64> vehicle_capacities{25,25,25,15,15,15};
      const std::vector<std::vector<int64>> groups{
                                           {0, 1, 2},//0
                                           {3, 4, 5},//1
                                           {6, 7, 8},//2
                                           {9,10,11}};//3
      const std::vector<RoutingIndexManager::NodeIndex > starts{
        RoutingIndexManager::NodeIndex{12},
        RoutingIndexManager::NodeIndex{12},
        RoutingIndexManager::NodeIndex{12},
        RoutingIndexManager::NodeIndex{12},
        RoutingIndexManager::NodeIndex{12},
        RoutingIndexManager::NodeIndex{12},
      };
      const std::vector<RoutingIndexManager::NodeIndex > ends{
        RoutingIndexManager::NodeIndex{13},
        RoutingIndexManager::NodeIndex{13},
        RoutingIndexManager::NodeIndex{13},
        RoutingIndexManager::NodeIndex{13},
        RoutingIndexManager::NodeIndex{13},
        RoutingIndexManager::NodeIndex{13},
      };
  }; //end_dada_model

  void PrintSolution(const DataModel& data, const RoutingIndexManager& manager,
                     const RoutingModel& routing, const Assignment& solution,
                     const RoutingDimension& time_dimension,
                     const RoutingDimension& capacity_dimension,
                     const std::vector<int64>  demand_service_duration) {
    int64 total_time{0};
    for (int vehicle_id = 0; vehicle_id < data.vehicle_capacities.size(); ++vehicle_id) {
      int64 index = routing.Start(vehicle_id);
      LOG(INFO) << "Route for vehicle " << vehicle_id << ":";
      std::ostringstream route;
      while (routing.IsEnd(index) == false) {
        auto time_var = time_dimension.CumulVar(index);
        auto capacity_var = capacity_dimension.CumulVar(index);
        auto idx = manager.IndexToNode(index).value();
        route << idx << " Time("
              << solution.Min(time_var) << ", " << solution.Max(time_var)
              << ")[" << demand_service_duration[idx]
              <<"]D(" << solution.Min(capacity_var) <<")-> ";
        index = solution.Value(routing.NextVar(index));
      }
      auto time_var = time_dimension.CumulVar(index);
      auto capacity_var = capacity_dimension.CumulVar(index);
      auto idx = manager.IndexToNode(index).value();
      LOG(INFO) << route.str() << idx << " Time("
                << solution.Min(time_var) << ", " << solution.Max(time_var)
                << ")[" << demand_service_duration[idx]
                <<"]D(" << solution.Min(capacity_var) <<")";
      LOG(INFO) << "Time of the route: " << solution.Min(time_var) << "min";
      total_time += solution.Min(time_var);
    }
    LOG(INFO) << "Total time of all routes: " << total_time << "min";
    LOG(INFO) << "";
    LOG(INFO) << "Advanced usage:";
    LOG(INFO) << "Problem solved in " << routing.solver()->wall_time() << "ms";
  }//end_pront_solution
//
  class BaseDimension {
  //base initial fot all dimensions
  public:
      RoutingDimension* base_dimension; // valiable for dimensions variables manipulations
      int64 base_transit_callback_idx; // store the dimension transit callback
      std::string name; //dimension name
      const RoutingIndexManager* manager; // the nodes manager
      RoutingModel* routing; // routing model reference structure
      BaseDimension(RoutingIndexManager* manager, RoutingModel* routing, std::string name)
      :
      manager(manager),routing(routing),name(name){}
  };

  class CapacityBaseDimension: public BaseDimension{
  public:
    std::vector<int64> capacity; //vehicles capacity
    std::vector<int64> demands;  // all demands

    CapacityBaseDimension(RoutingIndexManager* manager,
                          RoutingModel* routing,
                          std::vector<int64> capacity,
                          std::vector<int64> demands,
                          std::string name)
                          :
                          BaseDimension(manager, routing, name),
                          demands(demands), capacity(capacity){
      //Capacity callback
      RoutingTransitCallback2 demand_evaluator = [demands,manager](int64 from_index, int64 to_index) -> int64 {
                                                  // Convert from routingvariable Index to time matrix NodeIndex.
                                                  return demands[manager->IndexToNode(from_index).value()];};
      // add evaluator index
      base_transit_callback_idx = routing->RegisterTransitCallback(demand_evaluator);
      // add capacity dimension without slack
      routing->AddDimensionWithVehicleCapacity(base_transit_callback_idx,
                                              int64{0},
                                              capacity, //vehicles capacity
                                              true,
                                              name);
    }
  };

  class GroupContainsDimension: public BaseDimension{
  // control the group contains with callback function as io demand quantities
  public:
    const std::vector<std::vector<int64>> &io_quantities_matrix; 
    const std::vector<int64> quantity_horizon; // group capacity
    const std::vector<int64> demand_quantity;
    const int id;
    GroupContainsDimension(RoutingIndexManager* manager,
                              RoutingModel* routing,
                              const std::vector<std::vector<int64>> io_quantities_matrix,
                              const std::vector<int64> quantity_horizon,
                              std::string name,
                              const std::vector<int64> demand_quantity,
                              const int id)
                              :
                              BaseDimension(manager,routing, name),
                              io_quantities_matrix(io_quantities_matrix),
                              quantity_horizon(quantity_horizon),
                              demand_quantity(demand_quantity),
                              id(id){

      RoutingTransitCallback2 quantities_transit = [demand_quantity, io_quantities_matrix, manager]
                                                   (int64 from_index,int64 to_index)->int64{
                                                    // Convert from routing variable Index to io matrix NodeIndex.
                                                    auto i = manager->IndexToNode(from_index).value();
                                                    auto j = manager->IndexToNode(to_index).value();
                                                    // return the io direction (+1/-1) times demand quatity
                                                    return io_quantities_matrix[i][j]*std::abs(demand_quantity[j]);};
      // add evaluator index
      base_transit_callback_idx = routing->RegisterTransitCallback(quantities_transit);
      // add arc cost
      routing->SetArcCostEvaluatorOfAllVehicles(base_transit_callback_idx);
      // add base time dimension with slack
      routing->AddDimension(base_transit_callback_idx,
                            quantity_horizon[id],
                            quantity_horizon[id],
                            false,
                            name);
    }
  };

  class TimeBaseDimension: public BaseDimension{
  // this set include the time dimension with static transits
  public:
    const std::vector<std::vector<int64>> &time_matrix;
    const int64 time_horizon;
    const std::vector<int64> demand_time_duration;
    TimeBaseDimension(RoutingIndexManager* manager,
                      RoutingModel* routing,
                      const std::vector<std::vector<int64>> time_matrix,
                      const int64 time_horizon,
                      std::string nome,
                      const std::vector<int64> demand_time_duration) 
                      :
                      BaseDimension(manager,routing, nome),
                      time_matrix(time_matrix),
                      time_horizon(time_horizon),
                      demand_time_duration(demand_time_duration) {

      RoutingTransitCallback2 time_transit = [time_matrix, demand_time_duration, manager]
                                             (int64 from_index, int64 to_index) -> int64 {
                                             // Convert from routing variable Index to time matrix NodeIndex.
                                             auto i = manager->IndexToNode(from_index).value();
                                             auto j = manager->IndexToNode(to_index).value();
                                             // return the TT plus load/unload duration
                                             return time_matrix[i][j] + demand_time_duration[i];};
      // add evaluator index
      base_transit_callback_idx = routing->RegisterTransitCallback(time_transit);
      // add arc cost
      routing->SetArcCostEvaluatorOfAllVehicles(base_transit_callback_idx);
      // add base time dimension with slack
      routing->AddDimension(base_transit_callback_idx,
                            time_horizon,
                            time_horizon,
                            false,
                            name);
    }
  };
  
  class Spots{
  /*
    This set represents a group of nodes:
      - capacity  - is a variable for some quantity limits
      - contains  - is how much the group of nodes have of that quantity
      - locations - each id node in the spots
  */
  public:
    const int64 id, capacity, contains;
    const std::string name;
    const std::vector<int64>  locations;
    Spots(const int64 id, 
          const int64 capacity, 
          const int64 contains, 
          const std::string name, 
          const std::vector<int64>  locations)
          :
          id(id),
          capacity(capacity),
          contains(contains),
          name(name),
          locations(locations){}
  };

  class Terminal: public Spots{
  /*
    Each instance from this set represent some rules on spots.
    Govern in/out  parallel flow and provide queuing between
    vehicles for each demand. i.e., the constrains over nodes
    create vehicles behavior in order to model the terminal io.
    set_demand_visit() function use the or-tools variable IntervalVar
    to impose vehicle node visit
  */
  public:  
    const int64 number_of_par_lines; // {1} is and FIFO terminal, {2} the terminal have only two parallel lines ...
    std::vector<std::vector<int64>>   each_demand_time_per_line, // each vechicle time occuparion time
                                      each_demand_id_per_line;
    std::vector<IntervalVar*> intervals;
    const std::vector<std::pair<int64, int64>> v_serv;
    const std::vector<int64> vehicles; // all visiting vehicles
    Terminal(const int64 id, 
             const int64 capacity, 
             const int64 contains, 
             const std::string name, 
             const std::vector<int64>  locations,
             const std::vector<int64> vehicles,
             const int64 number_of_par_lines,
             const std::vector<std::pair<int64, int64>> v_serv)
             :
             Spots(id,capacity,contains,name,locations),
             vehicles(vehicles),
             number_of_par_lines(number_of_par_lines),
             v_serv(v_serv){}
    
    void configure_io(Solver* const solver,
                      const RoutingIndexManager manager,
                      RoutingModel* routing,
                      const RoutingDimension& time_dimension,
                      bool force){
      // this function impose the rules using contraints in
      // order to modeling the terminal-vehicles io interation
      //
      // setting the entry lines
      for (auto i = 0; i < this->number_of_par_lines; i++){
        this->each_demand_time_per_line.push_back({});
        this->each_demand_id_per_line.push_back({});
      }
      // input vehicles
      this->set_demand_visit(solver, manager, routing, time_dimension, force);
      // add the first demands to be met
      for(int i=0; i< this->each_demand_id_per_line.size();i++){
        this->each_demand_time_per_line[i].push_back(this->serv_fn(this->locations[i]));
        this->each_demand_id_per_line[i].push_back(i);
      }
      this->demand_occupation( this->locations.size() - 1, solver );
    }

    void demand_occupation(int demand,
                           Solver* const solver){
      // This build order of occupancy on terminal lines.
      // Add the fastest time occupation between lines, 
      // where the successor demand is imposed from the
      // first vehicle occupation time.
      // synchronization from the earliest time of the demands
      if  (demand > this->number_of_par_lines)
              demand_occupation( demand-1, solver);
      // the fastest queue
      auto min_line_id = this->find_min_time_between_lines();
      // add interval precedence constraint
      auto int1 = this->intervals[demand];
      auto int2 = this->intervals[
                    this->each_demand_id_per_line[min_line_id][each_demand_id_per_line[min_line_id].size() - 1]];
      solver->AddConstraint(
        solver->MakeEquality(this->intervals[demand]->StartExpr(),this->intervals[
                    this->each_demand_id_per_line[min_line_id][each_demand_id_per_line[min_line_id].size() - 1]]->EndExpr()));
      // mark the line with the interval demand alocation time
      this->each_demand_time_per_line[min_line_id].push_back(this->serv_fn(this->locations[demand]));
      // mark the line with id demand
      this->each_demand_id_per_line[min_line_id].push_back(demand);
    }

    int find_min_time_between_lines(){
      // The occupation for each parallel lines increases with
      //    "serv_fn", and this function return the shortest line time id.
      std::vector<int> min_line;
      int min_line_id = 0, min_val;
      for(int plines=0; plines < this->number_of_par_lines;plines++ ){
            min_val = std::accumulate(this->each_demand_time_per_line[plines].begin(),
                                      this->each_demand_time_per_line[plines].end(),
                                      0);
            min_line.push_back(min_val);
            if(min_val < min_line[plines-1])
                min_line_id = plines;
      }
      return min_line_id;
    }
 
    void set_demand_visit(Solver* const solver,
                          const RoutingIndexManager manager,
                          RoutingModel* routing,
                          const RoutingDimension& time_dimension,
                          bool force){
      //function use the or-tools variable IntervalVar
      //to impose vehicle node visit
      std::vector<int64> idx; 
      int demand = 0; // control the vehicle order
      for(int64 location:this->locations){
        idx.push_back(manager.NodeToIndex(RoutingIndexManager::NodeIndex(location)));
        if (force){
          this->force_vehicle_visit(solver, routing, idx[demand], this->vehicles[demand]);
        }else{
          //this->set_vehicle_visit();//TODO
        }
        this->intervals.push_back(
          solver->MakeFixedDurationIntervalVar(
            time_dimension.CumulVar(idx[demand]),
            this->serv_fn(location),
            "demand_interval"));
        demand++;
      }
    }

    void force_vehicle_visit(Solver* const solver, 
                             RoutingModel* routing,
                             int64 idx,
                             int64 vehicle){
      // Add constraint to force initialization node for a vehicle                               
      solver->AddConstraint(solver->MakeEquality(routing->VehicleVar(idx),routing->VehicleVar(routing->Start(vehicle))));
    }

    // void set_vehicle_visit(){}

    int64 serv_fn(int64 location){
      try
      {
        for (auto i=0; i < this->locations.size(); i++){
          if(this->locations[i] == this->v_serv[i].first){
            return this->v_serv[i].second;
          }
        }
      }
      catch(const std::exception& e)
      {
        std::cerr << e.what() << '\n';
      }
      return 0;
    }
  };

  void teste(){
    std::cout << "dentro do monitor " << std::endl;
  }

  void Run(){

    // Instantiate the data problem
    DataModel data;
    // Index manager
    RoutingIndexManager manager(data.time_matrix.size(),
                                data.vehicle_capacities.size(),
                                data.starts,
                                data.ends);
    // Routing model set
    RoutingModel routing(manager);
    // Cost for each vehicle
    // const std::vector<int64> cost;
    //
    // demand quantity times duration vector with pair set
    std::vector<int64>  demand_service_duration;
    for(auto i=0;i< data.demands.size();i++){
      demand_service_duration.push_back(0);
      if(data.demands[i] > 0){
        demand_service_duration[i] = data.demands[i]*data.demand_load_time[i]; //service total time
      }else{
        demand_service_duration[i] = data.demands[i]*data.demand_unload_time[i]*(-1); //service total time
      }
    }
    //
    // Add dimensions
    //
    // add capacity dimension for vehicle capacity solver control
    CapacityBaseDimension base_capacity_dimension(&manager,&routing,data.vehicle_capacities,data.demands,"VehiclesCapacity");
    const RoutingDimension& capacity_dimension = routing.GetDimensionOrDie(base_capacity_dimension.name);
    // add time-based transitions evaluations to support a dependent time dimension
    TimeBaseDimension base_time_dimension(&manager,&routing,data.time_matrix,data.time_horizon,"BaseTime",
                                          demand_service_duration);
    const RoutingDimension& time_dimension = routing.GetDimensionOrDie(base_time_dimension.name);
    // Add time window constraints for each location except depot.
    for (int i = 0; i < data.demands.size(); ++i) {
      int64 index = manager.NodeToIndex(RoutingIndexManager::NodeIndex(i));
      if(index != -1){        
        time_dimension.CumulVar(index)->SetRange(0,data.time_horizon);
      }
    }
    // Add time window constraints for each vehicle start node.
    for (int i = 0; i < data.vehicle_capacities.size(); ++i) {
      int64 index = routing.Start(i);
      time_dimension.CumulVar(index)->SetRange(0,data.time_horizon);
    }
    // add starts and ends routes for minimization
    for (int i = 0; i < data.vehicle_capacities.size(); ++i) {
      routing.AddVariableMinimizedByFinalizer(time_dimension.CumulVar(routing.Start(i)));
      routing.AddVariableMinimizedByFinalizer(time_dimension.CumulVar(routing.End(i)));
    }
    // add io quantities dimension
    std::vector<int64> terminal_capacity{10};
    // GroupContainsDimension io_quantities_dimension(&manager,&routing,
    //                                               data.io_demand_matrix,
    //                                               terminal_capacity,
    //                                               "IO_T0",
    //                                               data.demands,
    //                                               0); // id
    
    // const RoutingDimension& io_qnty_dimension = routing.GetDimensionOrDie(io_quantities_dimension.name);   
    // demand quantity times duration vector with pair set
    std::vector<std::pair<int64, int64>>  serv;
    for(auto i=0;i< data.groups[1].size();i++){
      serv.push_back(std::make_pair( // pair elements
                    data.groups[1][i], // location
                    data.demands[data.groups[1][i]]*data.demand_load_time[data.groups[1][i]])); //service total time
    }
    //
    //Add nodes groups restrictions
    //
    //Building a terminal for truck io
    Terminal PIP_truck = Terminal(0, //id
                                100, //capacity
                                 10, //contains
                        "PIP_truck", //name
                     data.groups[1], //locations
                            {3,4,5}, //vehicles ids
                                 3,  //number of parallel lines
                              serv); //vector with demand quantity times time duration 
    // set load rules                        
    PIP_truck.configure_io(routing.solver(),
                     manager,
                     &routing,
                     routing.GetDimensionOrDie("BaseTime"),
                     true); // impose visit                     
    std::vector<std::pair<int64, int64>>  serv_1;
    for(auto i=0;i< data.groups[3].size();i++){
      serv_1.push_back(std::make_pair( // pair elements
                    data.groups[3][i], // location
                    (-1)*data.demands[data.groups[3][i]]*data.demand_unload_time[data.groups[3][i]])); //service total time
    }
    //Set for trucks delivery
    Terminal PIP_truck_in = Terminal(1, //id
                                   100, //capacity
                                    10, //contains
                           "PIP_truck", //name
                        data.groups[3], //locations
                               {3,4,5}, //vehicles ids
                                     3,  //number of parallel lines
                                serv_1); //vector with demand quantity times time duration 
    // set in and unload rules                        
    PIP_truck_in.configure_io(routing.solver(),
                              manager,
                              &routing,
                              routing.GetDimensionOrDie("BaseTime"),
                              true); // impose visit                     
    //Building a terminal for trains io
    std::vector<std::pair<int64, int64>>  serv_2;
    for(auto i=0;i< data.groups[0].size();i++){
      serv_2.push_back(std::make_pair( // pair elements
                    data.groups[0][i], // location
                    data.demands[data.groups[0][i]]*data.demand_load_time[data.groups[0][i]])); //service total time
    }
    Terminal PIP_trains = Terminal(2, //id
                                 100, //capacity
                                  10, //contains
                        "PIP_trains", //name
                      data.groups[0], //locations
                             {0,1,2}, //vehicles ids
                                   1,  //number of parallel lines
                              serv_2); //vector with demand quantity times time duration 
    // set io rules                        
    PIP_trains.configure_io(routing.solver(),
                            manager,
                            &routing,
                            routing.GetDimensionOrDie("BaseTime"),
                            true); // impose visit
    //Building a terminal for trains io
    std::vector<std::pair<int64, int64>>  serv_3;
    for(auto i=0;i< data.groups[2].size();i++){
      serv_3.push_back(std::make_pair( // pair elements
                    data.groups[2][i], // location
                    data.demands[data.groups[2][i]]*data.demand_load_time[data.groups[2][i]])); //service total time
    }
    Terminal PIP_trains_in = Terminal(3, //id
                                    100, //capacity
                                     10, //contains
                           "PIP_trains", //name
                         data.groups[2], //locations
                                {0,1,2}, //vehicles ids
                                      1,  //number of parallel lines
                                 serv_3); //vector with demand quantity times time duration 
    // set io rules                        
    PIP_trains_in.configure_io(routing.solver(),
                     manager,
                     &routing,
                     routing.GetDimensionOrDie("BaseTime"),
                     true); // impose visit
    // Solver* const solver = routing.solver();
    // for (auto i=0; i<3; i++ ){
    //     solver->AddConstraint( //node precedence for trains transit validation
    //       solver->MakeEquality(
    //         PIP_trains.intervals[i]->EndExpr(), PIP_trains_in.intervals[i]->StartExpr()));
    // }
    Solver* const solver = routing.solver();
    LocalSearchMonitor* monitor = solver->GetLocalSearchMonitor();
    // DecisionBuilder* const db = solver->MakePhase();
    std::vector<SearchMonitor*> sm;
    
    
    // routing.AddSearchMonitor(monitor);
    // routing.AddAtSolutionCallback(teste);
    // add solver parameters
    RoutingSearchParameters searchParameters = DefaultRoutingSearchParameters();
    // set parameters
    searchParameters.set_first_solution_strategy(FirstSolutionStrategy::PATH_CHEAPEST_ARC);
    // get solution assignment
    const Assignment* solution = routing.SolveWithParameters(searchParameters);
    // print if exist a solution
    if (nullptr != solution){
      PrintSolution(data, manager, routing, *solution, time_dimension, capacity_dimension, demand_service_duration);
    }
    else{
      std::cout << "No solution!" << std::endl;
    }
   }
}


int main(){
     operations_research::Run();
#ifdef WIN32
	 system("pause");
#endif
     return EXIT_SUCCESS;
}