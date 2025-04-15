# NURBS

## Notes

- Everything we are required to do takes place in the `nurbs.cpp` file
- **Test the application** before submission! There is a IMGUI button in the app.
  - If the tests fail but the solution looks correctly, let the teachers know. There might be a precision error or something. These automated tests are kinda fuzzy.
- Don't put helper function into the application file, create them locally in the `nurbs` namespace I guess

## Hints

Everything in `nurbs.cpp`, function by function

### CORE FUNCTIONS

- `find_span` implement a binary search looking for (something?)
- `evaluate_basis_functions` from lecture part 1 - compute basis functions. p+1 basis functions is enough, no need for N+1
- `point_on_curve_in_homogeneous_space` this function already gets the computed basis functions. See slide 6 in presentation (section Hints)
- `point_on_surface_in_homogeneous_space` also slide 6.
- `point_on_surface_in_homogeneous_space` (the next one) copmute the knots in BOTH DIRECTIONS. Call the previous function to compute the point itself. This is the last function to implement to see the surface - after this the surface should be visible.
- `derivative_knot_vector` "very simple" - shortening the vector from sides. **This should be a one-liner**. Use the function hinted in the description (`assign`). 
- `derivative_control_grid_u` partial derivatives in the U direction. It operates on a vector inside of a vector. Operates on **rows**, it's easier to iterate. The vectors is a **list of rows**.
- `derivative_control_grid_v` partial derivatives in the V direction. Operates on **columns** - be careful how you access the data (the array of rows).
- `derivative_using_A` derivative of a rational surface, getting the final derivative vector. The `A` argument is basically `Cw` from slides - position in homogeneous space. The other one is a derivative - `Cw'`
- `derivative_of_surface_u` we now have all the data from the previous algorithms as arguments. This will be probably a one or two-liner of just calling the previously implemented functions with the correct data.
- `derivative_of_surface_v` same thing, but in a different direction (V)
- `derivative_of_surface_v` same thing, but in a different direction (V)

### Utility functions

Already implemented, don't worry about these.

## Troubleshoot

### Linux only - test cases

- There is a built-in bug. In the `test_cases/*.inc` files, in one of the first lines of the structure, there is
a float cast of a zero. That is supposed to be a boolean. Not sure if this applies also to the other `(float)0` appearances.

Currently it's like this:

```js
SurfaceGenerationParameters{ 5, 3, 4, 2, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)10, (float)0, false },
```

The line is supposed to be:

```js
SurfaceGenerationParameters{ 5, 3, 4, 2, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)10, (bool)0, false },
```

This causes problems only on linux. Windows compiler doesn't mind.
