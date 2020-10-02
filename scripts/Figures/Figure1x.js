//
// r2d3: https://rstudio.github.io/r2d3
//
const prior=data.filter(d=>d.type==="prior");
const posterior = data.filter(d=>d.type=="posterior");
//Margins

const margin = ({top: 20, right: 20, bottom: 50, left: 70});
 height = width/1.3;
// 

const yAxis = (text,scale)=>g => g
    .attr("transform", `translate(${margin.left-5},${margin.top})`)
    .attr("class","axis y")
    .call(d3.axisLeft(scale).ticks())
    .call(g => g.select(".tick:last-of-type text").clone()
          .attr("x", 20)
          .classed("axis-label",true)
          .attr("id","y-axis-label")
          .attr("text-anchor", "start")
          .attr("font-family", "sans-serif")
          .text(text))
          
const xAxis = (text,scale)=>g => g
    .attr("transform", `translate(${margin.left},${height-margin.bottom+5})`)
    .call(d3.axisBottom(scale).ticks(10))
    .call(g => g.append("text")
          .classed("axis-label",true)
          .attr("fill","currentColor")
          .attr("y", 40)
          .attr("x",((width-margin.right-margin.left)/2))
          .attr("text-anchor", "end")
          .text(text))



// Scales

const x=d3.scaleLinear().range([0,width-margin.left-margin.right])
          .domain(d3.extent(data,d=>d.Ne)).nice()
const y=d3.scaleLinear().range([height-margin.top-margin.bottom,0])
          .domain(d3.extent(data.map(d=>d.density)))
/// Make the figure 

  const group = svg.append("g")
      .attr("transform",`translate(${margin.left},${margin.top})`)
     
   svg.append("g")
      .call(xAxis("Ne",x));
 
  svg.append("g")
      .call(yAxis("density",y));
   
  svg.selectAll(".tick text").attr("font-size","18px")
  svg.selectAll(".axis-label").attr("font-size","18px")

  
  const p =  d3.line()
    .x(d => x(d.Ne))
    .y(d => y(d.density))
   
   const a = d3.area()
       .x(d => x(d.Ne))
       .y0(d => y(0))
       .y1(d=>y(d.density))
   
   group.append("path")
        .datum(posterior.filter(d=>d.fill!=="NA"))
          .attr("d",a)
          .attr("stroke-width",2)
          .attr("stroke","none")
          .attr("fill",d3.color(d3.schemeAccent[0]).brighter(0.2))
   
   group.append("path")
        .datum(posterior)
          .attr("d",p)
          .attr("stroke-width",3)
          .attr("stroke",d3.color(d3.schemeAccent[0]).darker(0.2))
          .attr("fill","none")
     group.append("path")
        .datum(prior)
          .attr("d",p)
          .attr("stroke-width",3)
          .attr("stroke",d3.color(d3.schemeAccent[1]).darker(0.2))
          .attr("fill","none")
   
    const scale = d3.scaleOrdinal()
  .domain(["Posterior", "Prior"])
  .range([d3.color(d3.schemeAccent[0]).darker(0.2),d3.color(d3.schemeAccent[1]).darker(0.2)]);
  


  
  
  svg.selectAll("text").attr("font-family","HelveticaNeue-Light, Helvetica Neue Light, Helvetica Neue, Helvetica, Arial, Lucida Grande, sans-serif") 