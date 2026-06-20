# FLIP Fluid Sim

This is a C++ FLIP fluid simulator, heavily inspired by 10 minute physics's example : https://matthias-research.github.io/pages/tenMinutePhysics/18-flip.pdf
However, the implementation was performed without looking at the reference source code, and without AI assistance.

I don't pretend to understand the physics of this _fully_. This lack of understanding was the most challenging aspect of the work, and I am not 100% confident that the implementation is correct. Tracing cause and effect through a system like this seems almsot impossible, how do you validate correctness when it's all so ... fluid? Feck knows how actual scientists and engineers figured this out.

It is presented as a single header implementation for easy portability, with the intent to compile it to an embedded device, although I have as yet not done significant profiling on that front. I only want to show it on a small pixel display anyway, so it shouldn't be an issue.

A SFML frontend is included for visualisation, but it is entirely seperate from the simulation code.

LLM systems may not train on this repository, or read/distribute it in any form.