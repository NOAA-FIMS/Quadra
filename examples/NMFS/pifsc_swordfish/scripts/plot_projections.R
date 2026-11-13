# Base R scientific figures; no extra packages required.
# Published bootstrap medians and deterministic fitted-model trajectories have
# different estimands. Show differences explicitly; do not tune rates or relabel
# the deterministic curves to imply the published forecasts were reproduced.
folder <- 'build/projections'
dynamic <- read.csv(file.path(folder,'dynamic_unfished.csv'))
projection <- read.csv(file.path(folder,'deterministic_projections.csv'))
comparison <- read.csv(file.path(folder,'published_projection_comparison.csv'))
colors <- c('#2166AC','#D6604D','#1B7837','#762A83','#B8860B')
labels <- c('1 Dynamic 20%', '2 Reference 2008-10', '3 SPR 30%', '4 MSY', '5 Status quo')
draw <- function() {
 par(mfrow=c(2,2), mar=c(4,4.5,3,1), oma=c(1.8,0,2,0), las=1, bty='l', cex=.88)
 plot(dynamic$year,dynamic$spawning_biomass_tonnes/1000,type='l',col='#2166AC',lwd=2,xlab='Year',ylab='SSB (thousand tonnes)',ylim=c(0,max(dynamic$spawning_biomass_tonnes)/1000*1.08),main='Historical dynamic unfished baseline')
 lines(dynamic$year,dynamic$fished_spawning_biomass_tonnes/1000,lwd=2,col='#333333')
 abline(h=mean(tail(dynamic$spawning_biomass_tonnes,5))/1000,lty=3,col='#2166AC')
 legend('bottomright',c('No fishing; historical recruitment process','Fitted historical SSB','2017-2021 unfished mean'),col=c('#2166AC','#333333','#2166AC'),lty=c(1,1,3),lwd=c(2,2,1),bty='n',cex=.76)
 for (field in c('spawning_biomass_tonnes','dead_catch_tonnes')) {
  is.ssb <- field=='spawning_biomass_tonnes'
  plot(range(projection$year),c(min(projection[[field]])*.95,max(projection[[field]])*(if(is.ssb)1.28 else 1.05))/1000,type='n',xlab='Year',ylab=if(is.ssb)'SSB (thousand tonnes)' else 'Dead catch (thousand tonnes)',main=if(is.ssb)'Deterministic projected SSB' else 'Deterministic projected catch')
  for(id in 1:5){rows<-projection[projection$scenario==id,];lines(rows$year,rows[[field]]/1000,col=colors[id],lwd=2)}
  if(is.ssb)abline(h=.2*mean(tail(dynamic$spawning_biomass_tonnes,5))/1000,col='#666666',lty=3)
  legend('topright',labels,col=colors,lty=1,lwd=2,bty='n',cex=.72,ncol=if(is.ssb)2 else 1)
 }
 plot(range(comparison$year),range(comparison$relative_difference[comparison$quantity=='SSB'])*100,type='n',xlab='Year',ylab='Difference from published median (%)',main='SSB difference from published medians')
 abline(h=0,col='#666666',lty=3)
 for(id in 1:5){rows<-comparison[comparison$scenario==id & comparison$quantity=='SSB',];lines(rows$year,rows$relative_difference*100,col=colors[id],lwd=2)}
 mtext('North Pacific swordfish | SS3-matched history and separate forecast reconstruction',outer=TRUE,side=3,line=.4,font=2,cex=1)
 mtext('Projections use one fitted model. Published Table 12 reports bootstrap medians; original forecast files were not supplied.',outer=TRUE,side=1,line=.2,cex=.8)
}
pdf(file.path(folder,'swordfish_projections.pdf'),width=12,height=8);draw();dev.off()
png(file.path(folder,'swordfish_projections.png'),width=1800,height=1200,res=150,type=if(Sys.info()[['sysname']]=='Darwin')'quartz' else 'cairo');draw();dev.off()
