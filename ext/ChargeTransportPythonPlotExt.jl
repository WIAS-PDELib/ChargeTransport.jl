module ChargeTransportPythonPlotExt
using PythonPlot
using ChargeTransport
using GridVisualize
function __init__()
    GridVisualize.default_plotter!(getproperty( ChargeTransportPythonPlotExt, :PythonPlot))
end
end
