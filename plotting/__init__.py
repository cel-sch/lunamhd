# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:25:41 2023

@author: celin
"""
from .growth_plotter import plot_growth
from .eigenfunc_plotter import plot_EF

Plotters = {}
Plotters["Growth"] = plot_growth
Plotters["EF"] = plot_EF
Plotters["Profiles"] = plot_profiles

__all__ = ["Plotters", "EF", "Profiles"]