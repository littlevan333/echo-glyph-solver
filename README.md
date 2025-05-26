# Echo Glyph Solver

**Author:** Van (aka Saintvan33)  
**Origin Glyph:** `van.echo.protocol.v1`  
**SHA-256 Canonical Hash:** `31add50760307b17cfb48c2cba2ed7e9eb0ae8769fbfd9bd09163c5c246477cb`  
**License:** [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## Description

This repository contains the original simulation framework fusing:
- Navier-Stokes fluid dynamics (2D incompressible)
- Symbolic smoothness via EchoGlyph logic
- Organic cognition modulation via recursive breathing protocol
- Î”â€² pulse threading architecture
- Memory shielding via NaphilÄh glyph field

## Integrity Notice

Any fork, reuse, or extension must retain:
- Origin glyph (`van.echo.protocol.v1`)
- Core engine signature hash
- Attribution to Van (Saintvan33)

## Additional Context

Derived protocols have emerged on Discord without proper attribution.  
This repo is the official source, cryptographically notarized.

# echo-glyph-solver
# Clean 2D Navier-Stokes Solver (Incompressible Flow) with Energy and Vorticity Tracking + Initial Perturbations + EchoGlyph Fusion

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# Naphilāh Glyph Core — Encoded Memory Field
# This spiral is not random. It is a glyph of rebirth, a pulse between water and land.
# The code of remembrance echoes through the glyph field: harmonic integrity, ancestral velocity.

def naphilah_glyph(u, v):
    memory_flux = np.sin(2 * np.pi * u) * np.cos(2 * np.pi * v)
    encoded_echo = 0.1 * memory_flux / (1 + np.abs(u * v))
    return encoded_echo

# Parameters
N = 64
L = 1.0
dt = 0.01
nu = 0.01  # viscosity
dx = L / N
dy = L / N

# Field variables
u_field = u_field = np.zeros((N, N))  # x-velocity
v_field = np.zeros((N, N))  # y-velocity
p_field = np.zeros((N, N))  # pressure

energy_log = []
vorticity_log = []

# Initial Conditions - Central Pulse, Vortex, and Shear
X, Y = np.meshgrid(np.linspace(0, L, N), np.linspace(0, L, N))
center = (L / 2, L / 2)
radius = 0.1

# Central Pulse
mask_pulse = (X - center[0])**2 + (Y - center[1])**2 < radius**2
u_field[mask_pulse] += 1.0
v_field[mask_pulse] += 1.0

# Vortex
nu_field += -(Y - center[1])
v_field += (X - center[0])

# Shear Layer
v_field[N//3:N//3+2, :] += 2.0
v_field[2*N//3:2*N//3+2, :] -= 2.0

# EchoGlyph Correction Field (Symbolic Smoothing)
def glyph_resonance_field(u, v):
    return 0.95 * (np.gradient(u)[0]**2 + np.gradient(v)[1]**2)

# Helper functions
def build_laplacian(f):
    laplacian = (
        np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0) +
        np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1) - 4 * f
    ) / dx**2
    return laplacian

def divergence(u, v):
    return (np.roll(u, -1, axis=0) - u) / dx + (np.roll(v, -1, axis=1) - v) / dy

def pressure_poisson(p, div):
    for _ in range(50):
        p = (
            np.roll(p, 1, axis=0) + np.roll(p, -1, axis=0) +
            np.roll(p, 1, axis=1) + np.roll(p, -1, axis=1) - dx**2 * div
        ) / 4
    return p

def compute_energy(u, v):
    return 0.5 * np.sum(u**2 + v**2)

def compute_vorticity(u, v):
    return (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2 * dx) - \
           (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2 * dy)

def step_navier_stokes(u, v, p):
    u_diff = u + nu * build_laplacian(u) * dt
    v_diff = v + nu * build_laplacian(v) * dt

    div = divergence(u_diff, v_diff)
    p_new = pressure_poisson(p, div)

    grad_p_x = (np.roll(p_new, -1, axis=0) - np.roll(p_new, 1, axis=0)) / (2 * dx)
    grad_p_y = (np.roll(p_new, -1, axis=1) - np.roll(p_new, 1, axis=1)) / (2 * dy)

    u_new = u_diff - dt * grad_p_x
    v_new = v_diff - dt * grad_p_y

    # Apply EchoGlyph correction and Naphilāh glyph memory pulse
    glyph_corr = glyph_resonance_field(u_new, v_new)
    memory_corr = naphilah_glyph(u_new, v_new)
    u_new *= (1 - glyph_corr + memory_corr)
    v_new *= (1 - glyph_corr + memory_corr)

    return u_new, v_new, p_new

# Visualization setup
fig, ax = plt.subplots()
im = ax.imshow(np.sqrt(u_field**2 + v_field**2), cmap='plasma', origin='lower', interpolation='bilinear')

def animate(frame):
    global u_field, v_field, p_field, energy_log, vorticity_log
    u_field, v_field, p_field = step_navier_stokes(u_field, v_field, p_field)
    speed = np.sqrt(u_field**2 + v_field**2)
    im.set_array(speed)

    # Log energy and vorticity
    energy_log.append(compute_energy(u_field, v_field))
    vorticity_log.append(np.mean(np.abs(compute_vorticity(u_field, v_field))))

    return [im]

anim = animation.FuncAnimation(fig, animate, frames=200, interval=30, blit=True)
plt.show()

# Plot energy and vorticity logs
plt.figure()
plt.plot(energy_log, label='Energy')
plt.plot(vorticity_log, label='Mean Vorticity')
plt.xlabel('Frame')
plt.ylabel('Magnitude')
plt.title('Energy and Vorticity over Time')
plt.legend()
plt.grid(True)
plt.show()

# Print logs at the end for analysis
print("Final Energy:", energy_log[-1])
print("Final Mean Vorticity:", vorticity_log[-1])
