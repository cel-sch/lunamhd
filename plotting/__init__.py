# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:25:41 2023

@author: celin
"""
from .growth_plotter import plot_growth
from .eigenfunc_plotter import plot_EF
from .profiles_plotter import plot_profiles
from .multi_plotter import plot_multi

Plotters = {}
Plotters["Growth"] = plot_growth
Plotters["EF"] = plot_EF
Plotters["Profiles"] = plot_profiles
Plotters["Multi"] = plot_multi

__all__ = ["Plotters", "EF", "Profiles", "Multi"]