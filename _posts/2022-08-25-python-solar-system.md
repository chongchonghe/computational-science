---
title: "Simulating the solar system with 70 lines of Python code"
author: "Chong-Chong He"
layout: post
comments: true
toc: false
image: images/cover-solar-system.jpg
hide: false
date: 2022-08-25
categories: [Python, animation]
---

Computational astrophysics is a really fun subject where we use computers to simulate astronomical objects and phenomena. To demonstrate the power and beauty of computer simulations in the study of astronomy, I show here a program I wrote with about 70 lines of Python code that simulates a (somewhat) realistic solar system.

We know that planets follow orbits determined mainly by the gravity from the Sun. To calculate the orbits of a planet, we need to calculate the gravitational force from the Sun and integrate it to get the velocity and position. We learned in general physics that the force following Newton's law of gravity combined with Newton's second law of motion gives the acceleration:
$$
\vec a = \frac{GM}{|\vec{r}|^2} \frac{\vec r}{|\vec{r}|}
$$
We use this formula to evolve the position and velocity of a planet:
$$
\vec{r}_{\rm new} = \vec{r}_{\rm old} + \vec{v}_{\rm old}~{\rm d}t\\
\vec{v}_{\rm new} = \vec{v}_{\rm old} + \vec{a}_{\rm new}~{\rm d}t
$$
This is a modified version of the Euler's method of integration, the backward Euler method. This method is stable and accurate enough for our purpose. 


The above equations can be written in component form:
$$
\begin{align}
\begin{bmatrix}
x \\ y \\ z
\end{bmatrix}
&= \begin{bmatrix}
x \\ y \\ z
\end{bmatrix}
+
\begin{bmatrix}
v_x \\ v_y \\ v_z
\end{bmatrix} \cdot {\rm d}t \\
\begin{bmatrix}
a_x \\ a_y \\ a_z
\end{bmatrix}
&=
\begin{bmatrix}
x \\ y \\ z
\end{bmatrix} \cdot \frac{GM}{(x^2+y^2+z^2)^{3/2}} \\
\begin{bmatrix}
v_x \\ v_y \\ v_z
\end{bmatrix}
&= \begin{bmatrix}
v_x \\ v_y \\ v_z
\end{bmatrix}
+ 
\begin{bmatrix}
a_x \\ a_y \\ a_z
\end{bmatrix} \cdot {\rm d}t
\end{align}
$$
The benefit of programming in Python is that Python handles vectors neatly. The nine equations above (three vector equations, each with three components) can be easily converted into just three lines of code in Python:

```python
p.r += p.v * dt
acc = -2.959e-4 * p.r / np.sum(p.r**2)**(3./2)  # in units of AU/day^2
p.v += acc * dt
```

Now we have the integrator, we need the initial conditions, i.e. the position and velocity of a planet at some time. Instead of making up initial conditions and creating a toy model, I choose to make use of the NASA JPL Horizons online solar system data to get the precise locations and velocities of a planet as an initial condition. (https://ssd.jpl.nasa.gov/?horizons) This is done by calling a Python module that sends queries to the Horizons database and obtains the data. This part of the code looks like this:

```python
from astroquery.jplhorizons import Horizons
obj = Horizons(id=1, location="@sun", epochs=Time("2017-01-01").jd, id_type='id').vectors()
```

where `id=1` corresponds to Mercury. Then, `obj['x']` gives the x position of Mercury on Jan 1st, 2017 (defined by "2017-01-01" in the second line).

Putting everything together, including animation via matplotlib animation tool, the complete code has only 75 lines of code. The program will generate the following animation. It computes the orbits of the 4 inner planets (Mercury, Venus, Earth, and Mars) in years 2019 and 2020. The code is available for download at this [GitHub repository](https://github.com/chongchonghe/Python-solar-system).

![solar_system_150dpi.mp4](images/solar_system_150dpi.mp4)

**What are real about this animation?**

The animation has accurate planet positions corresponding to the date shown at the top-left corner. The sizes of the objects are not to scale, although their relative sizes are correct.

**Does the code calculate the orbits in real-time, or it is just animating the data queried from the internet?**

It calculates the orbits in real-time. The initial conditions are from the internet. The program needs the positions and speeds of the planets at some date in order to predict their positions in future time.

**Why is the orbit of Mercury so ugly?**

Well, it is. Unlike any other planets in the solar system, Mercury has high orbital eccentricity, meaning its orbit is not very round.

**Why only 4 planets, not all 8 of them?**

Because starting from Jupiter, the 5th planet, the orbital radius becomes very large (>5 times that of Earth). Including them in the animation will make the inner planets nearly invisible.
