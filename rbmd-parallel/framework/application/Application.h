﻿#pragma once
#include "Configuration.h"
#include "JsonParser.h"
//#include "Object.h"
#include <memory>
#include <vector>

class CommandLine;
class System;
class Executioner;
class Output;
class ConsoleOutput;
class InitCondition;
class InitGlobal;
class HyperParameters;
class Execution;
class Para;

namespace vtkm
{
namespace cont
{
class DeviceAdapterId;
}
}

class Application
{
  //friend class MyTestFixture;

public:
  Application(int argc, char** argv);
  Application();
  virtual ~Application();

  virtual void PrintLogo();

  // 解析命令行
  virtual void ParseCLI();

  virtual void CreateActions();

  virtual void Run();
  //virtual void ExecuteCommandom() = 0;
  void SetupDevice();
  void OutputVersion();
  void HelpMessages();
  void ErrerMessages();

  void AddOutput(std::shared_ptr<Output> output);
  void AddInitCondition(std::shared_ptr<InitCondition> initCondition);

  auto& DeviceTag() { return this->_device; }
  auto& GetParser() { return this->_parser; }
  auto& GetSystem() { return this->_system; }
  auto& GetInitCondition() { return this->_init_condition; }
  auto& GetExecutioner() { return this->_executioner; }
  auto& GetOutputWarehouse() { return this->_owh; }
  auto& GetInitConditionWarehouse() { return this->_init_condition_wh; }
  auto& GetParameter() { return this->_parameter; }
  auto& GetRun() { return this->_run; }
  template<typename T>
  T& GetSystem()
  {
    static_assert(std::is_base_of<System, T>::value, "convert error");
    return static_cast<T&>(*_system);
  }

protected:
  //std::vector<std::shared_ptr<Action>> _awh; // action warehouse
  std::vector<std::shared_ptr<Output>> _owh; // output warehouse
  std::unique_ptr<CommandLine> _command_line;
  std::unique_ptr<InitGlobal> _init_global;
  std::shared_ptr<JsonParser> _parser;
  std::shared_ptr<System> _system;
  std::shared_ptr<Executioner> _executioner;
  std::shared_ptr<InitCondition> _init_condition;
  std::shared_ptr<Execution> _run;
  std::shared_ptr<Para> _parameter;
  //std::shared_ptr<Object> _obj;
  std::shared_ptr<Configuration> _cfg;
  std::shared_ptr<HyperParameters> _hyper_parameters;
  std::shared_ptr<Output> _Output;



  std::vector<std::shared_ptr<InitCondition>> _init_condition_wh; //init condition ware house




  // 不使用unique_ptr，因为需要#include DeviceAdapterId
  std::shared_ptr<vtkm::cont::DeviceAdapterId> _device;
};
