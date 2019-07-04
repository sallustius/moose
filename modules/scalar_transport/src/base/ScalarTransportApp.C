#include "ScalarTransportApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<ScalarTransportApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ScalarTransportApp::ScalarTransportApp(InputParameters parameters) : MooseApp(parameters)
{
  ScalarTransportApp::registerAll(_factory, _action_factory, _syntax);
}

ScalarTransportApp::~ScalarTransportApp() {}

void
ScalarTransportApp::registerAll(Factory & f, ActionFactory & af, Syntax & /*s*/)
{
  /* ModulesApp::registerAll(f, af, s); */
  Registry::registerObjectsTo(f, {"ScalarTransportApp"});
  Registry::registerActionsTo(af, {"ScalarTransportApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ScalarTransportApp::registerApps()
{
  registerApp(ScalarTransportApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ScalarTransportApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ScalarTransportApp::registerAll(f, af, s);
}
extern "C" void
ScalarTransportApp__registerApps()
{
  ScalarTransportApp::registerApps();
}
