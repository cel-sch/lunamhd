# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:25:41 2023

@author: celin
"""
from .growth_plotter import plot_growth

Plotters = {}
Plotters["Growth"] = plot_growth

__all__ = ["Plotters"]