# aslpoints

Simple experimental functions working on 2D and 3D point sets to fit planes, estimate geometric transforms, etc.

It is header-only currently. Requires the [ASL](https://github.com/aslze/asl) 1.11.10+ library.


Currently available functions (may change in the future):

* **fitPlaneXY**: fits a plane _z_ = _a_ _x_ + _b_ _y_ + _c_ to a set of 3D points and returns [a, b, c]

* **fitPlane**: fits a plane to a set of 3D points and returns it as a point and normal

* **fitPlane3D**: same with a different method (only one might survive)

* **fitCircle**: fits a circle to a set of 2D or 3D points

* **fitEllipse**: fits an ellipse to a set of 2D points

* **fitSphere**: fits a sphere to a set of 3D points.

* **findAffineTransform**: Estimates the affine transform between two sets of 2D or 3D points

* **findRigidTransform**: Estimates the rigid transform between two sets of 2D or 3D points

* **findHomography**: Computes a perspective transform between two sets of 2D points

* **fitPoly**: fits a polynomial (univariate y=f(x), bivariate z=f(x,y) or trivariate w=f(x,y,z)) of any degree to a set of 2D, 3D or 4D points and returns its coefficients

* **polynomial**: evaluates a polynomial (single variable, bivariate or trivariate)

Can easily be included using CMake in other projects where ASL is already. For example:



```cmake
include(FetchContent)
FetchContent_Declare(aslpoints URL https://github.com/aslze/aslpoints/archive/1.6.0.zip)
FetchContent_MakeAvailable(aslpoints)
target_link_libraries(myprogram asls aslpoints)
```
