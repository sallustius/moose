# Automatic Differentiation

## Objects Supported

!alert warning
Automatic differentiation is rapidly spreading through the framework, so it is
not a guarantee that this list will always be complete.

- `ADKernel`
- `ADVectorKernel`
- `ADNodalBC`
- `ADVectorNodalBC`
- `ADIntegratedBC`
- `ADVectorIntegratedBC`
- `ADDGKernel`
- `ADMaterial`

## Automatic Differentiation Migration Guide

### Top-level Declarations

* class and validParams predeclarations:

```c++
// from
class Foo

template <>
InputParameters validParams<Foo>();

// to
template <ComputeStage>
class Foo;

declareADValidParams(Foo);
```

* class declaration:

```c++
// from
class Foo : public Material/Kernel/BoundaryCondition/etc
{
  ...

// to
template <ComputeStage compute_stage>
class Foo : public ADMaterial/ADKernel/ADBoundaryCondition/etc<compute_stage>
{
  ...
```

* validParams implementation:

```c++
// from
template <>
InputParameters
validParams<Foo>()
{
  auto params = validParams<Material/Kernel/BoundaryCondition/etc>();
  params.addRequiredCoupledVar("bar", "blah");
  params.addParam("baz", "blah");
  ...
  return params;
}

// to
defineADValidParams(
    Foo,
    ADMaterial/ADKernel/ADBoundary/etc,

    params.addRequiredCoupledVar("bar", "blah");
    params.addParam("baz", "blah");
    ...
    );
```

* For AD intermediate base classes explicit instantiation is required if Foo inherits from an AD
 class and other classes inherit from Foo.  At the bottom of your `Foo.C` put:

```c++
adBaseClass(Foo);
```

### Materials

* variables (e.g. member vars):

```c++
// from
MaterialProperty<Foo>
// to
ADMaterialProperty(Foo)
```

* use material property:

```c++
// from
getMaterialProperty<Foo>(...);
// to
adGetADMaterialProperty<Foo>(...);
```

* declare material property:

```c++
// from
declareProperty<Foo>(...);
// to
adDeclareADProperty<Foo>(...);
```

### FE Variables and Coupling

* coupled stuff:

```c++
// from
coupledValue(...)
coupledGradient(...)
// to
adCoupledValue(...)
adCoupledGradient(...)
```

* Setting unused dimensions to zero

```c++
// something like this
for (unsigned i = _ndims; i < 3; ++i)
{
  _foo[i] = &adZero();
  _grad_foo[i] = &adGradZero();
}
```

### Variables

* predefined:

```
Real   -->   ADReal
RankTwoTensor   -->   ADRankTwoTensor
RankFourTensor   -->   ADRankFourTensor
???? more?
```
