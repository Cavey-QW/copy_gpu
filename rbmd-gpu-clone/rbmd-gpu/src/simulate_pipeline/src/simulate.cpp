#include "simulate.h"

int test_current_step = 0;
int test_num_steps;
Simulate::Simulate(std::shared_ptr<Ensemble>& simulate_pipeline,
                         std::shared_ptr<Output>& output)
    : _simulate_pipeline(simulate_pipeline)
    , _output(output)
{
    _time_step = DataManager::getInstance().getConfigData()->Get<rbmd::Real>("timestep", "execution");//0.001
    _num_steps = DataManager::getInstance().getConfigData()->Get<rbmd::Real>("num_steps", "execution");//10000
    test_num_steps = _num_steps;
    _current_time = 0;
    _current_step = 0;

}

void Simulate::Init() 
{    
    _simulate_pipeline->Init();
    _output->Init();
}

int Simulate::Execute() {
  while (KeepGoing()) {
    _output->Execute();
    _current_step++;
    _current_time += _time_step;
    test_current_step = _current_step;
    _simulate_pipeline->Run();
  }
  return 0;
}

bool Simulate::KeepGoing() {
  bool keep_going = false;

  if (_current_step <= _num_steps) keep_going = true;

  return keep_going;
}
