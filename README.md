# aslpoints

Simple experimental functions working on 2D and 3D point sets to fit planes, estimate geometric transforms, etc.

It is header-only currently. Requires the [ASL](https://github.com/aslze/asl) 1.11.8+ library.


Currently available functions (may change in the future):

* **fitPlaneXY**: fits a plane _z_ = _a_ _x_ + _b_ _y_ + _c_ to a set of 3D points and returns [a, b, c]

* **fitPlane**: fits a plane to a set of 3D points and returns it as a point and normal

* **findRigidTransform**: Estimates the rigid transform between two sets of 2D points

* **findRigidTransform**: Estimates the rigid transform between two sets of 3D points

* **findHomography**: Computes a perspective transform between two sets of 2D points

* **fitPoly**: fits a polynomial y=f(x) of any degree to a set of 2D points and returns its coefficients

* **fitPoly**: fits a bivariate polynomial z=f(x,y) to a set of 3D points and returns its coefficients



Can easily be included using CMake in other projects where ASL is already. For example:



```cmake
include(FetchContent)
FetchContent_Declare(aslpoints URL https://github.com/aslze/aslpoints/archive/1.0.zip)
FetchContent_MakeAvailable(aslpoints)
target_link_libraries(myprogram asls aslpoints)
```
