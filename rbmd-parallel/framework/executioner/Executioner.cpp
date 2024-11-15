﻿#include "Executioner.h"
#include "Application.h"
#include "Logging.h"
#include "Output.h"
#include "System.h"
#include <iostream>    
#include <iomanip>  

Executioner::Executioner(const Configuration& cfg)
  : Object(cfg)
  , _app(*(Get<Application*>("_app")))
  , _para(*(_app.GetParameter()))
  //, _system(*(_app.GetSystem()))
  , _run(*(_app.GetRun()))
  , _timer(_run.GetTimer())
  , _dt(_para.GetParameter<Real>(PARA_TIMESTEP))
  , _start_time(0)
  , _current_time(0)
  , _num_steps(_para.GetParameter<Real>(PARA_NUM_STEPS))
  , _end_time(std::numeric_limits<Real>::max())
{
  _current_step = 0;
  _end_time = _dt * _num_steps;
  _timer.Start();
}

void Executioner::Init() 
{
  //_system.Init();
  _run.Init();

  for (auto& output : _app.GetOutputWarehouse())
  {
    output->Init();
  }
}

void Executioner::Execute()
{
  PreExecute();

  while (KeepGoing())
  {
    Timer time;
    time.Start();
    //PreStep();
    TakeStep();
    PostStep();
  }

  PostExecute();
}

void Executioner::PostExecute()
{
  console::Success("运行完成.  运行时间:  ", _timer.GetElapsedTime());
  console::Info("当前步:", _current_step-1, "总步数:", _num_steps);
  console::Info("当前时间:", _current_time-_dt, "结束时间:",  _end_time);
  std::ofstream log_file("LogFile.csv", std::ios::app);
  try
  {
    log_file << "运行完成."<< std::endl <<"运行时间:   " << _timer.GetElapsedTime() << std::endl
             << "当前步:   " << _current_step-1 << std::endl << "总步数:   " << _num_steps << std::endl
             << "当前时间:   " << _current_time-_dt << std::endl<< "结束时间:   " << _end_time;
    log_file.close();
  }
  catch (const std::exception& e)
  {
    log_file.close();
    console::Error(e.what());
  }
}

void Executioner::PreStep() 
{
  ComputeDT();

  for (auto& output : _app.GetOutputWarehouse())
  {
    output->Execute();
  }
}

void Executioner::TakeStep() 
{
  //TODO: 求解前增加step，还是求解后增加step?
  _current_step++;
  std::cout << "当前步数：" << _current_step << std::endl;
  _current_time += _dt;
  _run.Execute();
  //_system.Evolve();
}

void Executioner::PostStep() {}

bool Executioner::KeepGoing()
{
  bool keep_going = false;

  if (_current_step <= _num_steps)
    keep_going = true;

  return keep_going;
}
