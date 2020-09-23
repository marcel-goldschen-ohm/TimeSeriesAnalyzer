# TimeSeriesAnalyzer
Flexible and performant viewer and analysis tool for groups of time series.

* Flexible grouping of multiple time series (e.g. recordings in multiple channels).
* Fast plotting of HUGE time series (requires Jim Hokanson's [plotBig_Matlab](https://github.com/JimHokanson/plotBig_Matlab)).
* Overlay each time series with any number of associated time series (e.g. idealization, fit, etc.)
    * A simple naming scheme allows the user to add an arbitrary number of these signals to the data.
* Dynamic selection of groups/data to visualize and traversal across time series within each group.
* Click and drag selection of x-axis range ROIs for analysis of subsections of time series.
    * !!! When active, all operations below will be restricted to selected subsections.
    * Functionality supplied by `XAxisROIManager` which can also be used standalone with any plot axes.
* Measure signal properties.
* Curve fitting.
* Mask, zero, interpolate, or apply basic mathematical operations.
* Simple idealization of piecewise continuous signals.
* Simple, flexible and easily extended underlying `struct array` data structure.
* Entire UI is in a single `uipanel` and is easily reparented into your own custom UI.

![User Interface](TimeSeriesAnalyzer.png "User Interface")
