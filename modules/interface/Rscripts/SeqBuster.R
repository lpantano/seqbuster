localhost<-"localhost"
user<-"lpantano"
db<-"seqbuster"
psswd<-"sqllorena"
require(tcltk)
# pathsource<-"http://172.18.11.210/seqbuster/"
pathsource<-paste(sep="",getwd(),"/")
#load database seqbuster
#source()
library(RMySQL)


MySQL(max.con = 16, fetch.default.rec = 5000000, force.reload = TRUE)

m <- dbDriver("MySQL")
# 
# con <- dbConnect(m,host=localhost,user=user,db=db,password=pssw)

showpackages<-function(){

# 	tkdelete(package)
# 	wpack<-tktoplevel()
# 	tkwm.geometry(wpack,"400x400+0+650")
	tkdelete(package,0,10)
	typeindex<-tkcurselection(type)
	typename<-tkget(type,typeindex)
	query<-paste(sep="","select `name` from `scriptsGUI` where `type` like '",typename,"';")
	print(query)
	rs <- dbSendQuery(con,query)
	packages<-as.data.frame(fetch(rs))
# 	print(packages)
	for (o in packages[,1]){
		tkinsert(package,"end",o)
	}
	
	

}
doanalysis<-function(){
	typeindex<-tkcurselection(package)
# 	print(typeidndex)
	typename<-tkget(package,typeindex)
# 	print(typename)
	packagename<-paste(collapse=" ",typename)
	query<-paste(sep="","select `script` from `scriptsGUI` where `name` like '",packagename,"';")
# 	print(query)
	rs <- dbSendQuery(con,query)
	script<-as.data.frame(fetch(rs))
	print(script[1,1])
	source(paste(sep="",pathsource,"GUIR/",script[1,1],"GUI.R"))
}


exitseqbuster<-function(){
	tkdestroy(wpro)
	tkdestroy(wtype)
	tkdestroy(wpack)
	tkdestroy(wsamples)
	dbDisconnect(con)
}

loadproject<-function(){
	
	pro<-tclvalue(project)
	query<-paste(sep="","select `name` from `project` where `name` like '",pro,"';")
	rs <- dbSendQuery(con,query)
	check<-as.data.frame(fetch(rs))
# 	print(project)
	if (nrow(check)>0){
		tkgrid(tklabel(wpro,text=paste("Load project:",pro)),column=1,row=1,columnspan=2)
# 		tkdelete(samples,0,10)
		query<-paste(sep="","select `name`,`description` from `",pro,"`.`experiments` ;")
# 	print(query)
		rs <- dbSendQuery(con,query)
		samplesdata<-as.data.frame(fetch(rs))
		# 		print(samplesdata)
		
		for (ind in 1:nrow(samplesdata)){
			
			tkgrid(tklabel(wsamples,text=paste(samplesdata[ind,1],":")),column=1,row=ind)
			tkgrid(tklabel(wsamples,text=samplesdata[ind,2]),column=2,row=ind)
		}
	
	}else{
		werror<-tktoplevel()
		tkwm.geometry(werror,"300x100+0+0")
		tkgrid(tklabel(werror,text="There is no project with that name"),column=1,row=2)
# 		tkgrid(tkbutton(werror,text="Load...",command=sub () tkdestroy(werror)),column=1,row=3)

	}
}

wpro<-tktoplevel()
tkwm.geometry(wpro,"200x100+0+0")
tkwm.title(wpro,"load project")
wtype<-tktoplevel()
tkwm.geometry(wtype,"300x300+0+200")
tkwm.title(wtype,"Analysis")
wpack<-tktoplevel()
tkwm.geometry(wpack,"300x300+0+550")
tkwm.title(wpack,"Packages")
wsamples<-tktoplevel()
tkwm.geometry(wsamples,"500x200+400+0")
tkwm.title(wsamples,"Samples")

project<-tclVar("testpro")
tkgrid(tklabel(wpro,text="Project:"),column=1,row=2)
tkgrid(tkentry(wpro,width=10,textvariable=project),column=2,row=2)

 con<<- dbConnect(m,host=localhost,user=user,db=db,password=psswd)
#con<<- dbConnect(m,host="localhost",user="lpantano",db="seqhand",password="sqllorena")
		
# 	tkdestroy(w)
	
	type<-tklistbox(wtype,height=6,selectmode="single",background="white",exportselection=0)
# 	options<-c("general","isomirs","clustering","differential expressino")
	query<-paste(sep="","select `type` from `scriptsGUI` group by `type`;")
	rs <- dbSendQuery(con,query)
	types<-as.data.frame(fetch(rs))

	for (o in types[,1]){
		tkinsert(type,"end",o)
	}
	tkgrid(tklabel(wtype,text="Select type of analysis"),column=1,row=2)
	tkgrid(type,column=1,row=3)
	tkgrid(tkbutton(wtype,text="Show packages",command=showpackages),column=1,row=4)
	
	package<-tklistbox(wpack,height=6,selectmode="single",background="white",exportselection=0)
	tkgrid(tklabel(wpack,text="Select package"),column=2,row=2)
	tkgrid(package,column=2,row=3)
	tkgrid(tkbutton(wpack,text="Do analysis",command=doanalysis),column=2,row=4)

tkgrid(tkbutton(wpro,text="Load...",command=loadproject),column=1,row=3)
tkgrid(tkbutton(wpro,text="Exit",command=exitseqbuster),column=2,row=3)


