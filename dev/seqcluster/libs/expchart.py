

HEAD='''

    <script src="../js/amcharts.js" type="text/javascript"></script>        
    <script type="text/javascript">
            var chart;

            %(data)s
function loadSeq (x){
    var chartData=allData[0][x];
    // SERIAL CHART
    chart = new AmCharts.AmSerialChart();
    chart.dataProvider = chartData;
    chart.categoryField = "sample";
    chart.startDuration = 1;
    
    // AXES
    // category
    var categoryAxis = chart.categoryAxis;
    categoryAxis.labelRotation = 45; // this line makes category values to be rotated
    categoryAxis.gridAlpha = 0;
    categoryAxis.fillAlpha = 1;
    categoryAxis.fillColor = "#FAFAFA";
    categoryAxis.gridPosition = "start";
    
    // value
    var valueAxis = new AmCharts.ValueAxis();
    valueAxis.dashLength = 5;
    valueAxis.title = "reads per sample";
    valueAxis.axisAlpha = 0;
    chart.addValueAxis(valueAxis);
    
    // GRAPH
    var graph = new AmCharts.AmGraph();
    graph.valueField = "reads";
    graph.colorField = "color";
    graph.balloonText = "<b>[[category]]: [[value]]</b>";
    graph.type = "column";
    graph.lineAlpha = 0;
    graph.fillAlphas = 1;
    graph.title = "abundance";
    chart.addGraph(graph);
    
    // CURSOR
    var chartCursor = new AmCharts.ChartCursor();
    chartCursor.cursorAlpha = 0;
    chartCursor.zoomable = false;
    chartCursor.categoryBalloonEnabled = false;
    chart.addChartCursor(chartCursor);
    
    
    // WRITE
    chart.write("chartdiv");
}



        </script>
'''

DIV='''
        <div id="chartdiv" style="width: 900px; height: 400px;"></div>
     <div class="footer">Created by armcharts.com<div>

'''

def addgraph(data):
    graph= HEAD % {'data':data}
    return graph


def getExpDiv():
    return DIV
#info=[{"DB":"Repeat","uni":34,"mul":45,"nocon":3},
#{"DB":"RNA","uni":15,"mul":4,"nocon":5}]
#bar=["uni","mul","nocon"]

#print createdata(info)
#print addgraph("unique","uni","green")
#print createchart(bar)




