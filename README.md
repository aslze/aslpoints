Simple experimental functions working on 2D and 3D point sets to fit planes and estimate geometric transforms.

It's a single header currently. Requires the ASL library. Can easily be included (FetchContent, submodule) in other projects where ASL is already.

```cmake
target_link_libraries(app asls aslpoints)
```
