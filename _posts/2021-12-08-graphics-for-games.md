---
title: Graphical simulation of an alien planet using Unity3D engine
date: 2021-12-08 15:46:22 +/-0000
categories: [Blogging, Demonstration]
tags: [unity, c#, hlsl, shaderlab]
---

## Overview
During this project, I built my skills in graphics programming by creating a graphical simulation of an alien planet surface.
The goal of the project was to integrate as many graphical effects as I can into the scene,
while keeping the frame-rate sufficiently high to maintain an interactive experience.
The specification was therefore deliberately open-ended.

![Window shadow](/assets/img/posts/graphics-for-games.gif){: .shadow style="max-width: 90%" }
_Cinematic scene images_

**Workflows used**: Git and Agile.
**Technologies used**: Unity, C#, HLSL, ShaderLab.

**Note**: *This project was done as part of assessment for module CSC3231.*

## Things Done
- A scene hierarchy where some elements move/rotate in relation to parent objects.
- Permanent change to the planet surface (meteor crashing into terrain).
- A navigable camera, moving around the 3D environment.
- Space-ship, clouds & northern lights graphical effects in the sky.
- A day/night cycle with appropriate lighting changes.
- A body of realistic looking water.
- Planet surface extends as far as the eye can see.
- Wind effect.

## YouTube Video of Project Running
<iframe width="560" height="315" src="https://www.youtube.com/embed/KdF_kYqG_yE" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Learning Outcomes
- Identified appropriate techniques for rendering graphics in real-time.
- Described graphical representations mathematically.
- Realised which advanced techniques are required to achieve realism.
- Used advanced techniques associated with lighting to create realism in graphical scenarios.
- Balanced processing and memory requirements of multiple graphical effects at cinematic frame rate.

## Learn More
For more knowledge about the project, visit my project's [repository](https://github.com/mantas2000/Graphics-For-Games){: target="_blank"}.
