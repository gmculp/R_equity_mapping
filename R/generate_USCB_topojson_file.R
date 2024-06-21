
###disable scientific notation###
options(scipen = 999)

###load packages###
library(data.table)
library(censusapi)
library(geos)
library(sf)
library(igraph)
library(geojsonio)

data.table::setDTthreads(1)

########################################################
###function using igraph to sort linestrings by nodes###
########################################################

###recursive graphing function###
r.ig_dt <- function(in.dt, t) {

	net <- graph_from_data_frame(in.dt, directed = TRUE)
	dg <- decompose(net)
	dg.dt <- rbindlist(lapply(1:length(dg),function(i){
		
		m.0 <- as_edgelist(dg[[i]])
		class(m.0) <- "numeric"
		dt.0 <- as.data.table(m.0)
		vec <- c()
		
		###extract first 3+ edge loop###
		vec <- as.numeric(names(girth(dg[[i]])$circle))
		
		###if loop exists###
		if(length(vec) > 0) {
			
			###girth function is not directed so check against original edges###
			dt.1a <- data.table(V1 = vec, V2 = c(vec[2:length(vec)],vec[1]))
			dt.1a[,o.id := .I]
			dt.1b <- data.table(V1 = c(vec[2:length(vec)],vec[1]), V2 = vec)
			dt.1b[,o.id := .I]
			setorder(dt.1b, -o.id)
			dt.1b[,o.id := .I]
			
			if(nrow(dt.1a[dt.0, on=.(V1,V2), nomatch = 0]) > nrow(dt.1b[dt.0, on=.(V1,V2), nomatch = 0])){
				dt.1 <- dt.1a[dt.0, on=.(V1,V2), nomatch = 0]
			} else{
				dt.1 <- dt.1b[dt.0, on=.(V1,V2), nomatch = 0]
			}
			
			rm(dt.1a,dt.1b)
			
			if(length(vec) > nrow(dt.1)){
				txt <- 'inc loop'
			} else{
				txt <- 'loop'
			}
			
			setorder(dt.1,o.id)
			dt.1[,o.id := NULL]
			dt.1[, sort.id := .I]
			dt.1[, sub.id := i+t]
			dt.1[, type := txt]
				
			###if loop is not entire graph###
			if(nrow(dt.1) < nrow(dt.0)){
			
				dt.2 <- r.ig_dt(dt.0[!dt.1, on=.(V1,V2)],t+1)
				this.out <- rbindlist(list(dt.1,dt.2),use.names=TRUE,fill=TRUE)
				#return(rbindlist(list(dt.1,dt.2),use.names=TRUE,fill=TRUE))

			###if loop is entire graph###	
			} else {
				this.out <- dt.1
				#return(dt.1)
			}
			
		###if graph is linear###	
		} else {
			xxx <- as.numeric(names(topo_sort(dg[[i]])))
			
			this.out <- data.table(V1=xxx[1:(length(xxx)-1)], V2=xxx[2:length(xxx)], sort.id=1:(length(xxx)-1), sub.id=i+t, type='line')

		}
		return(this.out)
	}),use.names=TRUE,fill=TRUE)
	return(dg.dt)
}


###function to generate graphs in parallel### 
sort_edges <- function(pairz, in.dt, f.node, t.node, in_clus=2){

	rel.dt <- data.table(pair.id=pairz)
	n.r <- nrow(rel.dt)	
	rel.dt[,pc := rep(1:in_clus,each=(ceiling(n.r/in_clus)))[1:n.r]]
	rel.dt[,r.id := seq_len(.N), by = list(pc)]

	out_dt <- data.table::rbindlist(parallel::mclapply(1:in_clus, function(zz) {
		temp.dt <- data.table::rbindlist(lapply(rel.dt[pc==zz,]$pair.id, function(j){
			t.dt <- r.ig_dt(in.dt[pair.id==j,c(f.node,t.node), with=FALSE], 0)
			t.dt[,pair.id := j]
			return(t.dt)
		}),use.names=TRUE,fill=TRUE)
		invisible(gc())
		return(temp.dt)
	}, mc.preschedule=FALSE, affinity.list = 1:in_clus),use.names=TRUE,fill=TRUE)
	
	setnames(out_dt, c('V1','V2'), c(f.node,t.node))
}




generate_USCB_topojson_file <- function(FIPS_dt, USCB_TIGER.path, geo.year="2020", geo.type="tract", omit.unpopulated=TRUE, omit.artifacts=TRUE, output.file_name, in_clus=2) {

	if(as.character(geo.year)=="2010") {
		f_year <- "2019"
	} else {
		f_year <- "2022"
		geo.year <- '2020'
	}
	
	FIPS.dt <- copy(as.data.table(FIPS_dt))
	
	###pad state and county codes with leading zeros###
	FIPS.dt[,state := sprintf("%02d", as.numeric(state))]
	FIPS.dt[,county := sprintf("%03d", as.numeric(county))]
	
	FIPS.dt <- unique(FIPS.dt[,c("state","county"),with=FALSE])
	
	########################
	###process face files###
	########################
	
	face.files <- paste0("tl_",f_year,"_",FIPS.dt$state,FIPS.dt$county,"_faces")
	
	###check if files are where they should be###
	ff <- list.files(file.path(USCB_TIGER.path,"FACES"))
	
	if(length(face.files[!face.files %in% ff]) > 0){
		stop("\nThe following files are missing:\n",paste(face.files[!face.files %in% ff],collapse="\n"))
	}

	faces.dt <- st_sf(rbindlist(lapply(face.files,function(j){

		###read in shapefile###
		temp.sf <- sf::st_read(file.path(USCB_TIGER.path,"FACES",j), j, stringsAsFactors = F, quiet=T)
		
		###convert to data.table###
		return(as.data.table(temp.sf))

	}), use.names=TRUE, fill=TRUE))

	faces.dt <- as.data.table(st_drop_geometry(faces.dt))

	###get ggeography columns###
	if(as.character(geo.year)=="2010") {
		faces.dt[,USCB_tract := paste0(STATEFP10,COUNTYFP10,TRACTCE10)]
		faces.dt[,USCB_block := paste0(STATEFP10,COUNTYFP10,TRACTCE10,BLOCKCE10)]
		faces.dt[,ZCTA := ZCTA5CE10]
		faces.dt[,PUMA := PUMACE10]
		faces.dt[,BG := BLKGRPCE10]
		
	} else {
		faces.dt[,USCB_tract := paste0(STATEFP20,COUNTYFP20,TRACTCE20)]
		faces.dt[,USCB_block := paste0(STATEFP20,COUNTYFP20,TRACTCE20,BLOCKCE20)]
		faces.dt[,ZCTA := ZCTA5CE20]
		faces.dt[,PUMA := PUMACE20]
		faces.dt[,BG := BLKGRPCE20]
	}
	
	########################
	###process edge files###
	########################

	edge.files <- paste0("tl_",f_year,"_",FIPS.dt$state,FIPS.dt$county,"_edges")
	
	###check if files are where they should be###
	ff <- list.files(file.path(USCB_TIGER.path,"EDGES"))
	
	if(length(edge.files[!edge.files %in% ff]) > 0){
		stop("\nThe following files are missing:\n",paste(edge.files[!edge.files %in% ff],collapse="\n"))
	}

	edges.sf <- st_sf(rbindlist(lapply(edge.files,function(j){

		###read in shapefile###
		temp.sf <- sf::st_read(file.path(USCB_TIGER.path,"EDGES",j), j, stringsAsFactors = F, quiet=T)
		
		###convert to data.table###
		return(as.data.table(temp.sf))

	}), use.names=TRUE, fill=TRUE))

	#TLID: Permanent edge ID
	#TFIDL: Permanent face ID on the left of the edge
	#TFIDR: Permanent face ID on the right of the edge
	#TNIDF: From TIGER node identifier
	#TNIDT: To TIGER node identifier


	###https://www.census.gov/programs-surveys/geography/about/glossary.html###
	###Block numbers beginning with a zero (in Block Group 0) are intended to include only water area###
	###but not all water-only blocks have block numbers beginning with 0 (zero).###

	temp.dt <- faces.dt[,c('TFID','USCB_block','BG','ZCTA','PUMA'),with=FALSE]

	temp.dt[,LWTYPE := ifelse(BG=='0','W','L')]
	temp.dt[,BG := NULL]
	temp.dt[,USCB_tract := substr(USCB_block,1,11)]

	setnames(temp.dt,names(temp.dt),paste0(names(temp.dt),'L'))
	edges.sf <- merge(edges.sf, temp.dt, by='TFIDL', all.x=TRUE)

	setnames(temp.dt,names(temp.dt),gsub('L$','R',names(temp.dt)))
	edges.sf <- merge(edges.sf, temp.dt, by='TFIDR', all.x=TRUE)
	rm(temp.dt)

	edges.sf$WKT <- geos_write_wkt(as_geos_geometry(st_geometry(edges.sf)))

	edges.sf$rem.flag <- ifelse(!is.na(edges.sf$USCB_blockL) & !is.na(edges.sf$USCB_blockR) & (edges.sf$USCB_blockL==edges.sf$USCB_blockR), 1, 0)

	###remove edges where both faces are water###
	edges.sf$rem.flag <- ifelse(!is.na(edges.sf$LWTYPEL) & !is.na(edges.sf$LWTYPER) & edges.sf$LWTYPEL=='W' & edges.sf$LWTYPER=='W', 1, edges.sf$rem.flag)

	edges.sf <- edges.sf[edges.sf$rem.flag==0,]
	edges.sf$rem.flag <- NULL

	###############################################
	###remove unpopulated contiguous land masses###
	###############################################

	mycensuskey <-"2ca0b2830ae4835905efab6c35f8cd2b3f570a8a"

	if(as.character(geo.year)=="2010") {
		my.survey <- "dec/sf1"
		my.vars <- c("P001001")
		my.vintage <- 2010
	} else {
		my.survey <- "dec/pl"
		my.vars <- c("P1_001N")
		my.vintage <- 2020
	}

	api.data_cb <- rbindlist(lapply(1:nrow(FIPS.dt), function(x) as.data.table(getCensus(name = my.survey,
		vintage = my.vintage,
		key = mycensuskey,
		vars = my.vars,
		region = "block:*",
		regionin = paste0("state:",FIPS.dt[x]$state,"+county:",FIPS.dt[x]$county)))))

	api.data_cb[,USCB_block := paste0(state,county,tract,block)]
	setnames(api.data_cb,my.vars,c("USCB_pop"))


	temp.dt <- unique(as.data.table(st_drop_geometry(edges.sf))[,c('USCB_blockL','LWTYPEL','USCB_blockR','LWTYPER'),with=FALSE])

	temp.dt[,USCB_blockR := ifelse(is.na(LWTYPER) | LWTYPER=='W',NA,USCB_blockR)]
	temp.dt[,USCB_blockL := ifelse(is.na(LWTYPEL) | LWTYPEL=='W',NA,USCB_blockL)]

	temp.dt[,USCB_block1 := pmin(USCB_blockR,USCB_blockL,na.rm=TRUE)]
	temp.dt[,USCB_block2 := pmax(USCB_blockR,USCB_blockL,na.rm=TRUE)]

	temp.dt <- unique(temp.dt[,c('USCB_block1','USCB_block2'),with=FALSE])

	net <- graph_from_data_frame(temp.dt[!is.na(USCB_block1) & !is.na(USCB_block2) & USCB_block1!=USCB_block2], directed = FALSE)
	dg <- decompose(net)
	dg.dt <- rbindlist(lapply(1:length(dg),function(i){
		data.table(USCB_block=as.character(names(V(dg[[i]]))),group=i)
	}))	

	dd <- unique(c(temp.dt$USCB_block1,temp.dt$USCB_block2))
	dd <- dd[!is.na(dd)]
	dd <- dd[!dd %in% dg.dt$USCB_block]

	dg.dt2 <- data.table(USCB_block=dd,group=((1:length(dd))+max(dg.dt$group)))


	unpop_contig.dt <- merge(rbindlist(list(dg.dt,dg.dt2),use.names=TRUE,fill=TRUE),api.data_cb[,c("USCB_block","USCB_pop"),with=FALSE],by="USCB_block",all.x=TRUE)

	unpop_contig.dt[,tot.pop:=sum(USCB_pop), by = group]

	unpop_contig.dt <- unpop_contig.dt[tot.pop==0]
	rm(dg,dg.dt,dg.dt2,net,dd,temp.dt,api.data_cb)


	###################################################################
	###determine artifacts misassigned to another county's waterbody###
	###################################################################

	temp.dt <- as.data.table(st_drop_geometry(edges.sf))[!is.na(USCB_blockL) & !is.na(USCB_blockR),.(tot=.N),by=list(USCB_blockL,LWTYPEL,USCB_blockR,LWTYPER)]

	temp.dt[,t.id := .I]

	###dedup R and L swaps###
	zing.dt <- merge(temp.dt,temp.dt,by.x=c('USCB_blockL','USCB_blockR'),by.y=c('USCB_blockR','USCB_blockL'))

	zing.dt[, t.id := pmin(t.id.x,t.id.y)]

	temp.dt <- temp.dt[!t.id %in% zing.dt$t.id]
	rm(zing.dt)

	temp.dt[,ctyL := substr(USCB_blockL,1,5)]
	temp.dt[,ctyR := substr(USCB_blockR,1,5)]

	art.dt <- rbindlist(list(temp.dt[LWTYPEL=='L' & !is.na(LWTYPER), c('USCB_blockL','LWTYPER','ctyR'), with=FALSE], temp.dt[LWTYPER=='L' & !is.na(LWTYPEL),c('USCB_blockR','LWTYPEL','ctyL'), with=FALSE]),use.names=FALSE)

	setnames(art.dt,c('USCB_blockL','LWTYPER','ctyR'),c('USCB_block','LWTYPE','cty'))

	art.dt <- art.dt[,.(tot=.N), by=.(USCB_block,LWTYPE,cty)]
	art.dt <- dcast(art.dt, USCB_block + cty ~ LWTYPE, value.var = "tot")

	art.dt <- merge(art.dt[substr(USCB_block,1,5) != cty & is.na(W) & !is.na(L)], art.dt[substr(USCB_block,1,5) == cty & !is.na(W) & is.na(L)], by="USCB_block")

	rm(temp.dt)



	#########################################################
	###limit to edges between polygons and coastline edges###
	#########################################################

	geom_edges.dt <- as.data.table(st_drop_geometry(edges.sf))

	if(grepl("P", geo.type, ignore.case=TRUE)){
		geo.type <- "PUMA"
		geom_edges.dt[,GEOID.R := PUMAR]
		geom_edges.dt[,GEOID.L := PUMAL]
	} else if(grepl("Z", geo.type, ignore.case=TRUE)){
		geo.type <- "ZCTA"
		geom_edges.dt[,GEOID.R := ZCTAR]
		geom_edges.dt[,GEOID.L := ZCTAL]
	} else {
		geo.type <- "USCB_tract"
		geom_edges.dt[,GEOID.R := USCB_tractR]
		geom_edges.dt[,GEOID.L := USCB_tractL]
	}
	

	###remove interior edges###
	geom_edges.dt[,rem.flag := ifelse(GEOID.R == GEOID.L & LWTYPER == LWTYPEL & !is.na(GEOID.R) & !is.na(GEOID.L) & !is.na(LWTYPER) & !is.na(LWTYPEL), 1, 0)]
	
	###remove unpopulated contiguous areas###
	if(omit.unpopulated){
		geom_edges.dt[,rem.flag := ifelse(!is.na(USCB_blockL) & USCB_blockL %in% unpop_contig.dt$USCB_block,1,rem.flag)]
		geom_edges.dt[,rem.flag := ifelse(!is.na(USCB_blockR) & USCB_blockR %in% unpop_contig.dt$USCB_block,1,rem.flag)]
	}

	###remove artifacts###
	if(omit.artifacts){
		geom_edges.dt[,rem.flag := ifelse(!is.na(USCB_blockR) & USCB_blockR %in% art.dt$USCB_block & !is.na(LWTYPEL) & LWTYPEL=='W',1,rem.flag)]
		geom_edges.dt[,rem.flag := ifelse(!is.na(USCB_blockL) & USCB_blockL %in% art.dt$USCB_block & !is.na(LWTYPER) & LWTYPER=='W',1,rem.flag)]
	}

	geom_edges.dt <- geom_edges.dt[rem.flag==0]
	geom_edges.dt[,rem.flag := NULL]

	geom_edges.dt[,GEOID.R := ifelse(!is.na(USCB_blockR) & USCB_blockR %in% art.dt$USCB_block,NA,GEOID.R)]
	geom_edges.dt[,GEOID.L := ifelse(!is.na(USCB_blockL) & USCB_blockL %in% art.dt$USCB_block,NA,GEOID.L)]

	###

	geom_edges.dt[,GEOID.R := ifelse(LWTYPER=='W',NA,GEOID.R)]
	geom_edges.dt[,GEOID.L := ifelse(LWTYPEL=='W',NA,GEOID.L)]

	geom_edges.dt <- unique(geom_edges.dt[,c("TLID","TNIDF","TNIDT","GEOID.L","GEOID.R","WKT"),with=FALSE])

	###remove NAs###
	geom_edges.dt <- geom_edges.dt[!(is.na(GEOID.L) & is.na(GEOID.R))]
	geom_edges.dt <- geom_edges.dt[!(is.na(TNIDF) & is.na(TNIDT))]

	###flip if is.na(GEOID.R)... occurrs at edge of study area###
	geom_edges.dt[,new.TNIDF := ifelse(is.na(GEOID.R),TNIDT,TNIDF)]
	geom_edges.dt[,new.TNIDT := ifelse(is.na(GEOID.R),TNIDF,TNIDT)]
	geom_edges.dt[,new.GEOID.L := ifelse(is.na(GEOID.R),GEOID.R,GEOID.L)]
	geom_edges.dt[,new.GEOID.R := ifelse(is.na(GEOID.R),GEOID.L,GEOID.R)]
	geom_edges.dt[,new.WKT := ifelse(is.na(GEOID.R), geos_write_wkt(geos_reverse(WKT)), WKT)]
	geom_edges.dt[,c("TNIDF","TNIDT","GEOID.L","GEOID.R","WKT") := NULL]

	setnames(geom_edges.dt,names(geom_edges.dt),gsub("^new\\.","",names(geom_edges.dt)))

	###sort by GEOID.L###
	setorder(geom_edges.dt, GEOID.L, GEOID.R)

	geom_edges.dt <- unique(geom_edges.dt)

	###assign unique ID starting at zero###
	geom_edges.dt[,u.id := .I] 

	geom_edges.dt[,u.id := u.id - 1]

	###############################################
	###generate graph for each geographic record###
	###############################################

	p1.dt <- geom_edges.dt[!is.na(GEOID.L), c('u.id','TNIDT','TNIDF','GEOID.L','WKT'), with=FALSE]
	p1.dt[,WKT := geos_write_wkt(geos_reverse(WKT))]

	###https://github.com/topojson/topojson-specification?tab=readme-ov-file#214-arc-indexes###
	###0 refers to the first arc, 1 refers to the second arc, and so on.
	###-1 refers to the reversed first arc, -2 refers to the reversed second arc, and so on
	p1.dt[,u.id := -(u.id+1)]


	p1.dt[,new.TNIDF := TNIDT]
	p1.dt[,new.TNIDT := TNIDF]
	p1.dt[,c("TNIDF","TNIDT") := NULL]
	setnames(p1.dt, c('GEOID.L',"new.TNIDF","new.TNIDT"), c('GEOID',"TNIDF","TNIDT"))

	p2.dt <- geom_edges.dt[!is.na(GEOID.R), c('u.id','TNIDT','TNIDF','GEOID.R','WKT'), with=FALSE]
	setnames(p2.dt, c('GEOID.R'), c('GEOID'))

	p.dt <- rbindlist(list(p1.dt,p2.dt),use.names=TRUE,fill=TRUE)
	rm(p1.dt,p2.dt)

	p.dt[, pair.id := .GRP, by=.(GEOID)]

	#
	##
	###


	###deal with self loops###
	e1.dt <- p.dt[TNIDF == TNIDT,c('u.id','pair.id','TNIDT','TNIDF'),with=FALSE]
	e1.dt[,type := "loop, one edge"]
	e1.dt[,sub.id := rowid(pair.id)]
	e1.dt[,sort.id := 1]
	setorder(e1.dt,pair.id,sub.id,sort.id)
	e1.dt[,c('TNIDT','TNIDF') := NULL]


	###deal with two edge loop###
	e2.dt <- merge(p.dt[TNIDF != TNIDT,c('u.id','pair.id','TNIDT','TNIDF'),with=FALSE],p.dt[TNIDF != TNIDT,c('u.id','pair.id','TNIDT','TNIDF'),with=FALSE],by.x=c('pair.id','TNIDT','TNIDF'),by.y=c('pair.id','TNIDF','TNIDT'))

	e2.dt <- e2.dt[u.id.x != u.id.y & u.id.x != -u.id.y]
	e2.dt[, r.id := .I]
	e2.dt <- melt(e2.dt, id = c("pair.id","r.id"), measure.vars = c("u.id.x","u.id.y"), value.name = "u.id")
	setorder(e2.dt,pair.id,r.id,-u.id)
	e2.dt[,variable := paste0('u.id.',rowid(r.id))]
	if(nrow(e2.dt) > 0){
		e2.dt <- dcast(e2.dt, pair.id + r.id ~ variable, value.var = "u.id")
		e2.dt[, r.id := NULL]
		e2.dt <- unique(e2.dt)
		e2.dt[, sub.id := rowid(pair.id)]
		e2.dt <- melt(e2.dt, id = c("pair.id","sub.id"), measure.vars = c("u.id.1","u.id.2"), value.name = "u.id")
		e2.dt[,sort.id := seq_len(.N), by = list(pair.id,sub.id)]
		e2.dt[,type := "loop, two edge"]
		setorder(e2.dt,pair.id,sub.id,sort.id)
		e2.dt[,variable := NULL]
		e3.dt <- p.dt[!(u.id %in% unique(c(e1.dt$u.id, e2.dt$u.id))) ,c('u.id','pair.id','TNIDT','TNIDF'),with=FALSE]
	} else{
		e3.dt <- p.dt[!(u.id %in% unique(e1.dt$u.id)) ,c('u.id','pair.id','TNIDT','TNIDF'),with=FALSE]
	}

	###deal with all other loops###
	
	
	e3.dt <- sort_edges(unique(e3.dt$pair.id), e3.dt, 'TNIDF', 'TNIDT', in_clus)
	e3.dt <- merge(e3.dt, p.dt[,c('u.id','pair.id','TNIDT','TNIDF'),with=FALSE], by=c('pair.id','TNIDF','TNIDT'))
	e3.dt[,c('TNIDT','TNIDF') := NULL]
	setorder(e3.dt,pair.id,sub.id,sort.id)

	###bind all loops together###
	e.dt <- rbindlist(list(e1.dt,e2.dt,e3.dt), use.names=TRUE, fill=TRUE)
	rm(e1.dt,e2.dt,e3.dt)
	setorder(e.dt,pair.id,type,sub.id,sort.id)

	###renumber sub.id so it is distict by pair.id###

	e.dt <- merge(e.dt, unique(e.dt[,c('pair.id','type','sub.id'), with=FALSE])[,new_sub.id := seq_len(.N), by = list(pair.id)], by=c('pair.id','type','sub.id'))

	setorder(e.dt, pair.id, type, sub.id, new_sub.id, sort.id)
	e.dt[,sub.id := new_sub.id]
	e.dt[,new_sub.id := NULL]

	p.dt <- merge(p.dt, e.dt, by=c('pair.id','u.id'), all.x=TRUE)
	setorder(p.dt,pair.id,sub.id,sort.id)
	rm(e.dt)

	#######################################################
	###create shapefile to check if geometries are valid### 
	#######################################################

	geom.dt <- p.dt[,.(WKT.line = geos_write_wkt(geos_line_merge(paste0("MULTILINESTRING (",paste(gsub("LINESTRING ","",WKT),collapse=", "),")"))), WKT.poly = geos_write_wkt(geos_unary_union(geos_polygonize(paste0("MULTILINESTRING (",paste(gsub("LINESTRING ","",WKT),collapse=", "),")")))), poly.arcs=paste0("[",paste(u.id,collapse=","),"]")), by=list(GEOID,pair.id,sub.id)]

	###check for invalid geometries###
	nn <- nrow(geom.dt[!grepl('^POLYGON',WKT.poly)])
	if(nn > 0){
		stop(paste("\nThere are",nn,"illegal geometries\n"))
	}
	
	###uncomment to view illegal geometries###
	#plot(st_geometry(st_as_sf(geos_read_wkt(geom.dt[!grepl('^POLYGON',WKT.line)]$out.WKT))))

	###assign unique feature id###
	geom.dt[,f.id := .I]

	###deal with holes###
	geom.dt[,cw := geos_is_clockwise(WKT.line)]


	###merge holes to polygons by pair.id and check if contained###
	holes.dt <- merge(geom.dt[!(cw)],geom.dt[(cw)],by='pair.id')
	holes.dt[,contains := geos_contains(WKT.poly.y,WKT.poly.x)]
	holes.dt <- holes.dt[(contains),.(tot.holes=.N, hole.WKT=paste(gsub("^POLYGON \\((.*)\\)$","\\1",WKT.poly.x),collapse=", "), hole.arcs=paste(poly.arcs.x,collapse=",")), by=.(f.id.y)]

	polys.dt <- merge(geom.dt[(cw)], holes.dt, by.x="f.id", by.y="f.id.y", all.x=TRUE)
	polys.dt[,out.WKT := ifelse(!is.na(hole.WKT), paste0(gsub("\\)$", "", WKT.poly), ",", hole.WKT,")"), WKT.poly)]
	polys.dt[,out.arcs := ifelse(!is.na(hole.WKT), paste0("[",poly.arcs,',',hole.arcs,"]"), paste0("[",poly.arcs,"]"))]

	###uncomment to preview###
	#plot(st_geometry(st_as_sf(geos_read_wkt(polys.dt$out.WKT))), col='gold')

	###create multipart polygons### 
	multi_poly.dt <- polys.dt[,.(WKT=paste0("MULTIPOLYGON (",paste(gsub("POLYGON ","",WKT.poly),collapse=", "),")"), arcs=paste0("[",paste(out.arcs,collapse=","),"]")), by = GEOID]

	setorder(multi_poly.dt,GEOID)
	multi_poly.dt[,poly.id := .I]

	###uncomment to preview###
	#plot(st_geometry(st_as_sf(geos_read_wkt(multi_poly.dt$WKT))), col='gold')

	########################
	###export as topojson###
	########################

	###collapse rows into arrays###
	multi_poly.dt[,str := paste0('{"type":"MultiPolygon","arcs":',arcs,',"properties":{"poly_id":',poly.id,',"GEOID":',GEOID,'}}')]

	###convert WKT string into nested array of coordinates###
	setorder(geom_edges.dt,u.id)
	geom_edges.dt[,arcs := gsub(", ","],[",WKT)]
	geom_edges.dt[,arcs := gsub("LINESTRING \\(","[",arcs)]
	geom_edges.dt[,arcs := gsub("\\)","]",arcs)]
	geom_edges.dt[,arcs := gsub(" ",",",arcs)]
	geom_edges.dt[,arcs := paste0("[",arcs,"]")]

	object_name <- paste0(geo.type,"_",geo.year)

	str.1 <- paste0('{"type":"Topology","objects":{"',object_name,'":{"type":"GeometryCollection","geometries":[')

	str.2 <- paste(multi_poly.dt$str,collapse=',')

	str.3 <- paste0(']}},"arcs":[',paste(geom_edges.dt$arcs,collapse=','),']')

	str.4 <- paste0(',"bbox":[',paste(as.numeric(wk::wk_bbox(geos_read_wkt(geom_edges.dt$WKT))),collapse=','),']}')

	fileConn <- file(output.file_name)
	writeLines(paste0(str.1,str.2,str.3,str.4), fileConn)
	close(fileConn)

	################################################
	###read in as sf and compare to multi_poly.dt###
	################################################
	
	###read in newly generated topojson file as sf object###
	topo.sf <- topojson_read(output.file_name)
	
	###uncomment to preview###
	#plot(st_geometry(topo.sf), col='gold')

	topo.sf$WKT <- geos_write_wkt(as_geos_geometry(st_geometry(topo.sf)))
	topo.dt <- as.data.table(st_drop_geometry(topo.sf))
	topo.dt[,GEOID := as.character(GEOID)]
	topo.dt <- merge(topo.dt, multi_poly.dt[,c('GEOID','WKT'), with=FALSE], by='GEOID', all.y=TRUE)
	
	fail.dt <- topo.dt[WKT.x != WKT.y]
	fail.dt[,is.valid := geos_is_valid(WKT.x)]
	
	nn <- nrow(fail.dt)
	n1 <- nrow(fail.dt[(is.valid)])
	n2 <- nrow(fail.dt[!(is.valid)])
	
	###uncomment to preview###
	#plot(st_geometry(st_as_sf(geos_read_wkt(topo.dt[WKT.x != WKT.y]$WKT.x))), col='gold')
	
	if(nn > 0){
		cat(paste("\nThere were geometry differences found in",nn,"polygons.\nOf these,",n1,"are valid and",n2,"are invalid.\n"))
	} else{
		cat("\nTopojson file successfully generated.\n")
	}
}
