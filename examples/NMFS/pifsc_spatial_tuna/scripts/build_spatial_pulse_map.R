args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else "build/assessment_outputs/data"
report_dir <- if (length(args) >= 2L) args[[2L]] else "build/assessment_outputs/report"

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("spatial pulse map requires the R package jsonlite")
}

spatial_path <- file.path(data_dir, "spatial_animation.csv")
if (!file.exists(spatial_path)) stop("missing spatial animation data: ", spatial_path)
spatial <- read.csv(spatial_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("record_type", "year", "season", "region", "fleet", "biomass",
              "retained_catch", "discard_catch", "from_region", "to_region",
              "movement_biomass", "movement_probability")
missing <- setdiff(required, names(spatial))
if (length(missing)) stop("spatial animation data missing columns: ", paste(missing, collapse = ", "))

projection_path <- file.path(data_dir, "projection_summary.csv")
projections <- if (file.exists(projection_path)) {
  read.csv(projection_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()

payload <- jsonlite::toJSON(
  list(spatial = spatial, projections = projections), dataframe = "rows",
  na = "null", auto_unbox = TRUE, digits = 12
)

html <- r"---(<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Tuna spatial fishery pulse</title>
<style>
:root{--ink:#dff7ff;--muted:#83aebd;--panel:#092333;--ocean:#041723;--cyan:#3ee6d0;--fleet1:#ffb84d;--fleet2:#f26d85}
*{box-sizing:border-box}body{margin:0;background:radial-gradient(circle at 50% 10%,#0c3a4b 0,#041721 52%,#020c12 100%);color:var(--ink);font:15px system-ui,-apple-system,sans-serif;min-height:100vh}
main{max-width:1220px;margin:auto;padding:28px}.head{display:flex;justify-content:space-between;gap:24px;align-items:end}.eyebrow{text-transform:uppercase;letter-spacing:.18em;color:var(--cyan);font-size:11px;font-weight:800}h1{font-size:clamp(28px,4vw,48px);margin:5px 0 3px;line-height:1}.sub{color:var(--muted)}
.card{background:linear-gradient(145deg,rgba(9,35,51,.94),rgba(4,23,35,.88));border:1px solid #1a4a5d;border-radius:18px;box-shadow:0 18px 60px #0008;margin-top:20px;overflow:hidden}.toolbar{display:flex;align-items:center;gap:10px;flex-wrap:wrap;padding:15px 18px;border-bottom:1px solid #164052}.toolbar button{border:1px solid #2b6578;background:#0a3142;color:var(--ink);border-radius:999px;padding:8px 14px;cursor:pointer;font-weight:700}.toolbar button.active,.toolbar button:hover{background:#12647a}.toolbar input[type=range]{flex:1;min-width:180px;accent-color:var(--cyan)}.stamp{font-variant-numeric:tabular-nums;font-weight:800;min-width:130px;text-align:right}
#ocean{display:block;width:100%;height:auto;background:radial-gradient(ellipse at center,#0d4353 0,#062532 65%,#041923 100%)}.route{fill:none;stroke:#43e6d2;stroke-linecap:round;opacity:.48;filter:drop-shadow(0 0 6px #38d9c4)}.flowdot{fill:#d8fff8;filter:drop-shadow(0 0 5px #fff)}.region{fill:#0d5667;stroke:#75fff0;stroke-width:2;filter:drop-shadow(0 0 18px #21bca7)}.region-label{fill:#fff;font-size:19px;font-weight:800;text-anchor:middle}.metric{fill:#a8d1dc;font-size:12px;text-anchor:middle}.fleet-ring{fill:none;stroke-linecap:round;opacity:.92}.legend{display:flex;gap:24px;flex-wrap:wrap;padding:14px 20px;color:var(--muted)}.swatch{width:10px;height:10px;border-radius:50%;display:inline-block;margin-right:7px}
.race{padding:18px 20px 22px;border-top:1px solid #164052}.race h2{margin:0 0 2px;font-size:19px}.race p{margin:0 0 10px;color:var(--muted)}#race{width:100%;height:auto;display:block}.axis{stroke:#477181;stroke-width:1}.scenario{fill:none;stroke-width:4;stroke-linecap:round;stroke-linejoin:round}.scenario-label{font-size:12px;font-weight:800}.empty{fill:#83aebd;font-size:16px;text-anchor:middle}.foot{color:#6f98a6;font-size:12px;margin-top:12px}
@media(max-width:650px){main{padding:14px}.head{display:block}.stamp{text-align:left}.toolbar input[type=range]{order:9;flex-basis:100%}}
</style></head><body><main>
<div class="head"><div><div class="eyebrow">Quadra assessment explorer</div><h1>Spatial fishery pulse</h1><div class="sub">Biomass, seasonal movement, fleet catch, and management futures</div></div><div class="stamp" id="headline"></div></div>
<section class="card"><div class="toolbar"><button id="play">▶ Play</button><button id="history" class="active">Historical</button><button id="future">Scenario race</button><label>Speed <select id="speed"><option value="1200">Slow</option><option value="650" selected>Normal</option><option value="260">Fast</option></select></label><input id="time" type="range" min="0" value="0"><span class="stamp" id="stamp"></span></div>
<svg id="ocean" viewBox="0 0 1000 480" role="img" aria-label="Animated regional biomass and fishing activity"></svg>
<div class="legend" id="legend"><span><i class="swatch" style="background:#3ee6d0"></i>movement biomass</span><span>Node size = regional biomass</span><span>Ring width = fleet catch; solid = retained, faint = discard</span></div>
<div class="race"><h2>Management strategy race</h2><p>Projected spawning biomass by fishing scenario. Select “Scenario race” to animate the futures.</p><svg id="race" viewBox="0 0 1000 300"></svg></div></section>
<div class="foot">Generated from spatial_animation.csv and projection_summary.csv. Regions are schematic because the synthetic assessment does not assign geographic coordinates.</div>
</main><script>
const DATA=__DATA__;
const S=DATA.spatial||[], P=DATA.projections||[];
const regionRows=S.filter(d=>d.record_type==="region"), moveRows=S.filter(d=>d.record_type==="movement");
const frames=[...new Set(regionRows.map(d=>`${d.year}-${d.season}`))].map(k=>{const [year,season]=k.split("-").map(Number);return{year,season}}).sort((a,b)=>a.year-b.year||a.season-b.season);
const regions=[...new Set(regionRows.map(d=>+d.region))].sort((a,b)=>a-b), fleets=[...new Set(regionRows.map(d=>+d.fleet))].sort((a,b)=>a-b);
const fleetColors=["#ffb84d","#f26d85","#8ca8ff","#c69cff","#f5ee72"];
const maxBio=Math.max(1,...regionRows.map(d=>+d.biomass||0)), maxCatch=Math.max(1,...regionRows.map(d=>(+d.retained_catch||0)+(+d.discard_catch||0)));
let mode="history", index=0, timer=null;
const ocean=document.querySelector("#ocean"), slider=document.querySelector("#time"), stamp=document.querySelector("#stamp"), headline=document.querySelector("#headline");
const esc=s=>String(s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[c]));
const fmt=x=>Number(x||0).toLocaleString(undefined,{maximumFractionDigits:0});
function positions(){return Object.fromEntries(regions.map((r,i)=>[r,{x:180+(640*(i+.5)/regions.length),y:235+(i%2?35:-35)}]));}
function historicalFrame(i){const f=frames[i], rows=regionRows.filter(d=>+d.year===f.year&&+d.season===f.season), moves=moveRows.filter(d=>+d.year===f.year&&+d.season===f.season), pos=positions();let out="";
 moves.filter(d=>+d.from_region!==+d.to_region).forEach((d,j)=>{const a=pos[+d.from_region],b=pos[+d.to_region],w=2+18*Math.sqrt((+d.movement_biomass||0)/maxBio),bend=(+d.from_region<+d.to_region?-1:1)*(55+18*j);const path=`M${a.x},${a.y} Q${(a.x+b.x)/2},${(a.y+b.y)/2+bend} ${b.x},${b.y}`;out+=`<path id="route${j}" class="route" stroke-width="${w}" d="${path}"/><circle class="flowdot" r="4"><animateMotion dur="1.8s" repeatCount="indefinite"><mpath href="#route${j}"/></animateMotion></circle>`});
 regions.forEach(r=>{const p=pos[r],rr=rows.filter(d=>+d.region===r),bio=+(rr[0]?.biomass||0),radius=42+50*Math.sqrt(bio/maxBio);out+=`<circle class="region" cx="${p.x}" cy="${p.y}" r="${radius}"/>`;
 rr.forEach((d,j)=>{const retained=+d.retained_catch||0,discard=+d.discard_catch||0,ring=radius+10+j*10,width=2+12*Math.sqrt((retained+discard)/maxCatch),color=fleetColors[j%fleetColors.length];out+=`<circle class="fleet-ring" cx="${p.x}" cy="${p.y}" r="${ring}" stroke="${color}" stroke-width="${width}" stroke-dasharray="${retained+discard?`${Math.max(8,100*retained/(retained+discard))} ${Math.max(3,100*discard/(retained+discard))}`:"1 99"}"/>`});
 const totalCatch=rr.reduce((v,d)=>v+(+d.retained_catch||0)+(+d.discard_catch||0),0);out+=`<text class="region-label" x="${p.x}" y="${p.y-5}">Region ${r}</text><text class="metric" x="${p.x}" y="${p.y+17}">${fmt(bio)} biomass</text><text class="metric" x="${p.x}" y="${p.y+34}">${fmt(totalCatch)} catch</text>`});
 const total=regions.reduce((v,r)=>v+(+(rows.find(d=>+d.region===r)?.biomass)||0),0);headline.textContent=`Total biomass ${fmt(total)}`;stamp.textContent=`Year ${f.year} · Season ${f.season}`;ocean.innerHTML=out;}
const scenarios=[...new Set(P.map(d=>d.scenario))], years=[...new Set(P.map(d=>+d.projection_year))].sort((a,b)=>a-b), scenarioColors=["#3ee6d0","#ffb84d","#f26d85","#8ca8ff","#c69cff"];
function raceFrame(limit=years.length-1){const svg=document.querySelector("#race");if(!P.length||!years.length){svg.innerHTML='<text class="empty" x="500" y="145">Projection results are not available for this run.</text>';return}const vals=P.map(d=>+d.spawning_biomass||0),max=Math.max(...vals,1),x=y=>70+850*(years.indexOf(y)/Math.max(1,years.length-1)),yy=v=>250-205*v/max;let out='<line class="axis" x1="70" y1="250" x2="940" y2="250"/><line class="axis" x1="70" y1="35" x2="70" y2="250"/>';
 scenarios.forEach((s,i)=>{const rows=P.filter(d=>d.scenario===s&&years.indexOf(+d.projection_year)<=limit).sort((a,b)=>+a.projection_year-+b.projection_year);if(!rows.length)return;const points=rows.map(d=>`${x(+d.projection_year)},${yy(+d.spawning_biomass)}`).join(" "),last=rows[rows.length-1],color=scenarioColors[i%scenarioColors.length];out+=`<polyline class="scenario" stroke="${color}" points="${points}"/><circle fill="${color}" cx="${x(+last.projection_year)}" cy="${yy(+last.spawning_biomass)}" r="6"/><text class="scenario-label" fill="${color}" x="${Math.min(945,x(+last.projection_year)+9)}" y="${yy(+last.spawning_biomass)-7}">${esc(s)}</text>`});
 out+=`<text class="metric" x="500" y="286">Projection year</text><text class="metric" transform="translate(18 150) rotate(-90)">Spawning biomass</text>`;svg.innerHTML=out;}
function render(){if(mode==="history"){slider.max=Math.max(0,frames.length-1);index=Math.min(index,frames.length-1);historicalFrame(index);raceFrame(years.length-1)}else{slider.max=Math.max(0,years.length-1);index=Math.min(index,years.length-1);raceFrame(index);stamp.textContent=`Projection year ${years[index]||1}`;headline.textContent="Management futures"}}
function stop(){clearInterval(timer);timer=null;document.querySelector("#play").textContent="▶ Play"}function play(){if(timer){stop();return}document.querySelector("#play").textContent="❚❚ Pause";timer=setInterval(()=>{const max=+slider.max;index=index>=max?0:index+1;slider.value=index;render()},+document.querySelector("#speed").value)}
document.querySelector("#play").onclick=play;slider.oninput=e=>{index=+e.target.value;render()};document.querySelector("#speed").onchange=()=>{if(timer){stop();play()}};
function setMode(m){stop();mode=m;index=0;slider.value=0;document.querySelector("#history").classList.toggle("active",m==="history");document.querySelector("#future").classList.toggle("active",m==="future");render()}
document.querySelector("#history").onclick=()=>setMode("history");document.querySelector("#future").onclick=()=>setMode("future");
fleets.forEach((f,i)=>document.querySelector("#legend").insertAdjacentHTML("beforeend",`<span><i class="swatch" style="background:${fleetColors[i%fleetColors.length]}"></i>Fleet ${f}</span>`));render();
</script></body></html>)---"

parts <- strsplit(html, "__DATA__", fixed = TRUE)[[1L]]
html <- paste0(parts[[1L]], payload, parts[[2L]])
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(report_dir, "spatial_fishery_pulse.html")
writeLines(html, output_path, useBytes = TRUE)
cat("Spatial fishery pulse map written to", output_path, "\n")
