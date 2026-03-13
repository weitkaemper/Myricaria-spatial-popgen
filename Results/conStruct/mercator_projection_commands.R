deg_to_rad <- function(x){return (x * 2 * pi/360)}
projection <- function(y){return (log(tan(0.25*pi +0.5*deg_to_rad(y))))}


library(conStruct)
load(_results_)
load(_data block_)


mercator <- cbind(deg_to_rad(data.block$coords[,1]), projection(data.block$coords[,2]))
colnames(mercator) <- colnames(data.block$coords)

pdf("mercator_DN2_chain5_pies.pdf")
make.admix.pie.plot(admix.proportions = conStruct.results$chain_5$MAP$admix.proportions, coords = mercator,  x.lim = c(min(mercator[,1] - deg_to_rad(1)),max(mercator[,1] + deg_to_rad(1))), y.lim = c(min(mercator[,2] - deg_to_rad(1)),max(mercator[,2] + deg_to_rad(1))))
dev.off()
