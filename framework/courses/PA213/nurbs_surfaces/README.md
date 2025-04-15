# NURBS


## Troubleshoot

- There is a built-in bug. In the `test_cases/*.inc` files, in one of the first lines of the structure, there is
a float cast of a zero. That is supposed to be a boolean.

Currently it's like this:

```js
SurfaceGenerationParameters{ 5, 3, 4, 2, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)10, (float)0, false },
```

The line is supposed to be:

```js
SurfaceGenerationParameters{ 5, 3, 4, 2, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)10, (bool)0, false },
```

This causes problems only on linux. Windows compiler doesn't mind.
