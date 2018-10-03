/* License: GNU General Public License
SeqBuster: uder-friendly tool for small RNA analysis. Copyright (C) 2009 Lorena Pantano
 This program is free software: you can redistribute it and/or modify it under the terms
 of the GNU General Public License as published by the Free Software Foundation, either
 version 3 of the License, or (at your option) any later version. This program is distributed
 in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 more details. You should have received a copy of the GNU General Public License along with this
 program. If not, see http://www.gnu.org/licenses/.
 */


package miraligner;
import java.io.*;
import java.util.Date;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 *
 * @author Lorena Pantano
 */

public class map {

    private static Logger LOGGER = Logger.getLogger(Logger.GLOBAL_LOGGER_NAME);
    public static void readseq (String namein,String namedb,String sp,int mm,int tri,int add,String f,String nameo, boolean freq, boolean precursor, Integer minl) throws FileNotFoundException, IOException{
        LOGGER.info(new Date()+"\n");
        
        String l="";
        int annotate=0;
        Integer namecode=0;
        String last="";
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(nameo+".mirna")));
        PrintWriter outp = new PrintWriter(new BufferedWriter(new FileWriter(nameo+".mirna.opt")));
        PrintWriter outnm = new PrintWriter(new BufferedWriter(new FileWriter(nameo+".mirna.nomap")));
        TreeMap<String,TreeMap<Integer,Integer>> clusters = new TreeMap<String,TreeMap<Integer,Integer>>();
        TreeMap<String,TreeMap<String,Integer[]>> micropos = new TreeMap<String,TreeMap<String,Integer[]>>();
        TreeMap<Integer,Integer> tempclus = new TreeMap<Integer,Integer>();
        TreeMap<Integer,String> hashseq = new TreeMap<Integer,String>();
        TreeMap<Integer,String> nameseq = new TreeMap<Integer,String>();
        TreeMap<String,String> listmirna =  new TreeMap<String,String>();
        TreeMap<String,String> listinfo =  new TreeMap<String,String>();

        outp.printf("Initiated at: %s\n",new Date());
        outp.printf("Mismatches allowed: %s\n",mm);
        outp.printf("Addition allowed: %s\n",add);
        outp.printf("Trimming allowed: %s\n",tri);
        outp.close();
        TreeMap<Integer,TreeMap<String,alignment>> mapinfo = new TreeMap<Integer,TreeMap<String,alignment>>();
        TreeMap<Integer,Double> scoreseq = new TreeMap<Integer,Double>();
        TreeMap<String,String> preseq = new TreeMap<String,String>();
        
        BufferedReader inmi= new BufferedReader(new FileReader(namedb+"/miRNA.str"));
         while ((l=inmi.readLine())!=null){
            if (l.contains(">")){
                String[] pre=l.split(" ");
                pre[0]=pre[0].replace(">", "");
                l=l.replaceAll("\\]",",");
                l=l.replaceAll("\\[",",");
                String [] name=l.split(",");
                for (int i=1;i<name.length;i++){
                  if (name[i].contains("-")){
                    String [] namepos=name[i].split(":");
                    String [] pos=namepos[1].split("-");
                    Integer [] posi= new Integer[2];
                    posi[0]= Integer.parseInt(pos[0]);
                    posi[1]= Integer.parseInt(pos[1]);
                    if (micropos.containsKey(pre[0])){
                        TreeMap<String,Integer[]> mi=micropos.get(pre[0]);
                        mi.put(namepos[0],posi);
                        micropos.put(pre[0],mi);
                    }else{
                        TreeMap<String,Integer[]> mi= new TreeMap<String,Integer[]>();
                        mi.put(namepos[0],posi);
                        micropos.put(pre[0],mi);
                    }
                    
                  }
                }
            }
        }

        inmi.close();
       
        TreeMap<String,String>  reads=tools.getseq(namein,f);
        
        for (String name : reads.keySet()) {
            l=reads.get(name);
            namecode++;
            nameseq.put(namecode,name);
            if (l.length()>=minl & l.length()<=50){
                l=l.replaceAll("U","T");
                hashseq.put(namecode,l);
                String seed1=l.substring(0,8);
                String seed2=l.substring(8,16);
                
                if (clusters.containsKey(seed1)){

                    tempclus=clusters.get(seed1);
                    tempclus.put(namecode, 1);
                    clusters.put(seed1,tempclus);
                }else{

                    TreeMap<Integer,Integer> newclus = new TreeMap<Integer,Integer>();
                    newclus.put(namecode,1);
                    clusters.put(seed1,newclus);
                }
                if (clusters.containsKey(seed2)){

                    tempclus=clusters.get(seed2);
                    tempclus.put(namecode, 2);
                    clusters.put(seed2,tempclus);
                }else{

                    TreeMap<Integer,Integer> newclus = new TreeMap<Integer,Integer>();
                    newclus.put(namecode,2);
                    clusters.put(seed2,newclus);
                }
            }
        }
        
        LOGGER.info("Number of reads to be mapped: "+namecode);

        int bi=0;
        int ei=4,numlines=0;
        last="";
        String name="";
        int pospre=0;
        
        LOGGER.info("Searching in precursors");
        BufferedReader indb= new BufferedReader(new FileReader(namedb+"/hairpin.fa"));
        while ((l=indb.readLine())!=null){
            if (!l.contains(">") & name.contains(sp)){
                l=l.replaceAll("U","T");
                if (preseq.containsKey(name)){

                    String seqt=preseq.get(name);
                    seqt=seqt+l;
                    preseq.put(name, seqt);
                }else{
                     String seqt=l;
                     preseq.put(name, seqt);
                }
            }else{

                String [] pre=l.split(" ");
                name=pre[0];
                name=name.replaceAll(">","");
                
            }

        }
        for (String namechr : preseq.keySet()) {
                l=preseq.get(namechr);
                l=l.toUpperCase();
                pospre=1;
                numlines++;
                for (int i=bi;i<l.length()-8;i++){
                    ei=i+8;
                    String nt=l.substring(i, ei);
                    if (clusters.containsKey(nt)){
                        tempclus=clusters.get(nt);
                        for (int codeseq : tempclus.keySet()) {
                            //go to alignments
                            int pospretemp=pospre+i;
                            String seqdb=l.substring(i, l.length())+"NNN";
                            int posseed=tempclus.get(codeseq);
                            String seq=hashseq.get(codeseq);
 
                            if ((posseed==1 ) | (posseed==2 & i>=8)){
                                if (posseed==2 ){
                                    seqdb=l.substring(i-8, l.length())+"NNN";
                                    pospretemp=pospre+i-8;
                                }
                                if (seqdb.length()>=seq.length()){
                                alignment alg=align2(hashseq.get(codeseq),seqdb,posseed);
                                LOGGER.finest("->mapped: "+namechr+" position:"+pospretemp+" score: "+alg.scmut+" changes: "+alg.mut+" addition: "+alg.add);
                                if (alg.scmut<=mm & alg.add.length()<=add){
                                    LOGGER.finest("-->It is good.");
                                    alg.pospre=pospretemp;
                                    //do micro alignment
                                    if (!scoreseq.containsKey(codeseq)){
                                        scoreseq.put(codeseq,alg.sc);
                                        TreeMap<String,alignment> prealg = new TreeMap<String,alignment>();
                                        alg.amb=1;
                                        prealg.put(namechr,alg);
                                        mapinfo.put(codeseq, prealg);
                                    }else{
                                        if (alg.sc==scoreseq.get(codeseq)){
                                            TreeMap<String,alignment> prealg = mapinfo.get(codeseq);
                                            if (!prealg.containsKey(namechr)){
                                             alg.amb++;
                                            }
                                            prealg.put(namechr,alg);
                                            mapinfo.put(codeseq, prealg);
                                        }else if (alg.sc>scoreseq.get(codeseq)){
                                            scoreseq.put(codeseq,alg.sc);
                                            mapinfo.remove(codeseq);
                                            TreeMap<String,alignment> prealg = new TreeMap<String,alignment>();
                                            alg.amb=1;
                                            prealg.put(namechr,alg);
                                            mapinfo.put(codeseq, prealg);
                                        }
                                    }
                                }
                                }
                            }
                        }
                    }
                }

        }


        annotate=mapinfo.size();
        indb.close();
        if (freq){
            out.printf("seq\tname\tfreq\tmir\tstart\tend\tmism\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\tambiguity\n");
        }else{
            out.printf("seq\tname\tmir\tstart\tend\tmism\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\tambiguity\n");
        }
        LOGGER.finest("filtering output");
        for (int nc : mapinfo.keySet()) {
            TreeMap<String,alignment> pret=mapinfo.get(nc);
            int overlapp=0;
            String seq=hashseq.get(nc);
            LOGGER.finest("->name: "+nameseq.get(nc)+" sequence: " +seq);
            int ambmir=pret.size();
            listmirna.clear();
            listinfo.clear();
            for (String p : pret.keySet()) {
                overlapp=0;
                alignment at=pret.get(p);
                int end=at.pospre+at.hit-1;
                LOGGER.finest("-->precursor: "+p+" start: "+at.pospre+" end: "+end+" matches: "+at.hit+" changes: "+at.mut);
                if (micropos.containsKey(p)){
                TreeMap<String,Integer[]> mi = micropos.get(p);
                

                for (String m : mi.keySet()) {
                   int overlap3=0;
                   int overlap5=0;
                   Integer [] pos=mi.get(m);
                   LOGGER.finest("--->reference: "+m+" start: "+pos[0]+" end: "+pos[1]);
                   //overlap
                   if (pos[0]-at.pospre<=tri & pos[0]-at.pospre>0){
                        //get t5
                        at.t5=preseq.get(p).substring(at.pospre-1,pos[0]-1).toUpperCase();
                        LOGGER.finest("---->q5 template addition: "+at.t5);
                        overlap5=1;
                   }
                   if (at.pospre-pos[0]<=tri & at.pospre-pos[0]>0){
                        //get t5
                        at.t5=preseq.get(p).substring(pos[0]-1,at.pospre-1).toLowerCase();
                        LOGGER.finest("---->t5 truncation: "+at.t5);
                        overlap5=1;
                   }

                   if (pos[1]-end<=tri & pos[1]-end>0){
                        //get t3
                        at.t3=preseq.get(p).substring(end,pos[1]).toLowerCase();;
                        LOGGER.finest("---->t3 truncation: "+at.t3);
                        overlap3=1;
                   }
                   if (preseq.get(p).length()>end){
                    if (end-pos[1]<=tri & end-pos[1]>0 ){
                        //get t3
                        at.t3=preseq.get(p).substring(pos[1],end).toUpperCase();;
                        LOGGER.finest("---->q3 template addition: "+at.t3);
                        overlap3=1;
                    }
                   }else{
                       at.t3="longer";
                   }
                   if ( at.pospre-pos[0]==0){
                        
                        overlap5=1;
                        at.t5="0";
                   }
                   if(end-pos[1]==0){
                        overlap3=1;
                        at.t3="0";
                   }
                   LOGGER.finest("--->is inside miRNA range: "+overlap3+" "+overlap5);
                   if (overlap3==1 & overlap5==1 ){
                       overlapp=1;
                       int min=5;
                       int max=4;
                       if (pos[0]-5<0){min=0;}
                       if (pos[1]+4>preseq.get(p).length()){max=preseq.get(p).length()-pos[1];}
                       at.s5=preseq.get(p).substring(pos[0]-min,pos[0]+3);
                       at.s3=preseq.get(p).substring(pos[1]-4,pos[1]+max);
                       String ann="";
                       if (freq){
                           ann=seq+"\t"+nameseq.get(nc)+"\t"+tools.getFreq(nameseq.get(nc))+"\t"+m+"\t"+at.pospre+"\t"+end+"\t"+at.mut+"\t"+at.add.replace("u-", "")+"\t"+at.t5+"\t"+at.t3+"\t"+at.s5+"\t"+at.s3+"\tmiRNA\t"+p+"\t";      
                       }else{
                           ann=seq+"\t"+nameseq.get(nc)+"\t"+m+"\t"+at.pospre+"\t"+end+"\t"+at.mut+"\t"+at.add.replace("u-", "")+"\t"+at.t5+"\t"+at.t3+"\t"+at.s5+"\t"+at.s3+"\tmiRNA\t"+p+"\t";
                       }      
                       
                       listinfo.put(seq+m,ann);
                       listmirna.put(m,seq);
                        
                       LOGGER.finest("--->to mirna: " + seq+"\t"+nameseq.get(nc)+"\t"+m+"\t"+at.pospre+"\t"+end+"\t"+at.mut+"\t"+at.add+"\t"+at.t5+"\t"+at.t3+"\t"+at.s5+"\t"+at.s3+"\tmiRNA");
                       
                     }
                }
            }
            if (overlapp==0 & precursor==true){
               String pann="";
               if (freq){
                   pann=hashseq.get(nc)+"\t"+nameseq.get(nc)+"\t"+tools.getFreq(nameseq.get(nc))+"\t"+p+"\t"+at.pospre+"\t"+end+"\t"+at.mut+"\t"+at.add.replace("u-", "") + "\t0\t0\t0\t0\tprecursor\tNA\t"+ambmir+"\n";      
               }else{
                   pann=hashseq.get(nc)+"\t"+nameseq.get(nc)+"\t"+p+"\t"+at.pospre+"\t"+end+"\t"+at.mut+"\t"+at.add.replace("u-", "") +"\t0\t0\t0\t0\tprecursor\tNA\t"+ambmir+"\n";
               }
                  out.printf(pann);
                  LOGGER.finest("--->to precursor: "+seq+"\t"+nameseq.get(nc)+"\t"+p+"\t"+at.pospre+"\t"+end+"\t"+at.mut+"\t"+at.add+"\tNA\tNA\tNA\tNA\tprecursor\t"+ambmir);      
            }
            }
            
            ambmir=listmirna.size();
            for (String mirann: listinfo.keySet()){
                out.printf(listinfo.get(mirann)+ambmir+"\n");
            }

        }
        
        for (int nc:hashseq.keySet()){
            if (!mapinfo.containsKey(nc)){
                outnm.printf(">seq%s\n%s\n",nc,hashseq.get(nc));

            }

        }

        out.close();
        outnm.close();
        clusters.clear();
        tempclus.clear();
        mapinfo.clear();
        hashseq.clear();
        micropos.clear();
        scoreseq.clear();
        preseq.clear();

        LOGGER.info("" + new Date());
        LOGGER.info("Num reads annotated: "+annotate);
    }


    public static int substitution (String a,String b,int penalty,int reward){
        int score=penalty;
        if (a.matches(b)){
            score=reward;
        }
        return score;

    }
    private static alignment align2(String seq, String db,int seed){
        alignment alg=new alignment();
        int sc=0,mism=0,lastscore=0;
        int ns=0,lp=0,seedmism=0;
        for (int p=0;p<=seq.length()-11;p+=8){

            ns++;
            lp=p;
            String seed2=seq.substring(p,p+8);
            String seeddb2=db.substring(p,p+8);
            if (seed2.matches(seeddb2)){
                sc++;
            }else{
                seedmism=p;
            }
        }
        lp+=8;
        if (ns-sc==1){
            //go to look for mism
            String seed2=seq.substring(seedmism,seedmism+8);
            String seeddb2=db.substring(seedmism,seedmism+8);
            alg=alignment(seed2,seeddb2,alg,seedmism);
            if (alg.scmut==1){
                String t=seq.substring(lp,seq.length());
                String tdb=db.substring(lp,lp+t.length());
                alg=align3end(t,tdb,mism,alg,lp);
            }else{
               alg.scmut=10;
            }
        }else if (ns-sc==0){
            //go to trimming
            String t=seq.substring(lp,seq.length());
            String tdb=db.substring(lp,lp+t.length());
            alg=align3end(t,tdb,mism,alg,lp);
        }else{
            alg.scmut=10;
        }
        alg.hit=seq.length()-alg.add.length();
        if (alg.add.equals("0")){
            alg.hit += 1;
        }
        alg.sc=seq.length()-alg.scmut-alg.add.length()*0.5;
        return alg;
    }


    public static alignment align3end (String seq,String db,int mm,alignment alg,int pos){
        String [] seqNT=seq.split("");
        String [] adNT=db.split("");
        int score=0;
        int minlen=seq.length()-3;
        alg=alignment(seq.substring(0,minlen),db.substring(0,minlen),alg,pos);
        for (int i=seq.length()-1;i>minlen-1;i--){
            int sc=substitution(seqNT[i],adNT[i],0,1);
            if (sc==0){
                mm++;
                score=i+1;
            }
        }
        if (mm>0){
            alg.add=seq.substring(score-1,seq.length());
        }
        return alg;
    }
    public static alignment alignment (String seq,String db,alignment alg,int pos){
        int score=0;
        int i=0;

        int g;
        g=-3;
        int mm=0;
        int pe=0;
        int re=1;
        String [] seqNT=seq.split("");
        String [] adNT=db.split("");
        int minlen;
        minlen=seq.length();
        for (i=0;i<minlen;i++){
            int sc=substitution(seqNT[i],adNT[i],pe,re);
            if (sc==0){
                alg.scmut++;
                int tpos=i+pos+1;
                alg.mut=tpos+seqNT[i]+""+adNT[i];
            }
            score+=sc;
            if (i>score+2){
                score=0;
                i=minlen+10;
            }
        }
        return alg;
     }

    public static boolean removefiles (String str,String pathout){

    File dir = new File(pathout);
    String[] children = dir.list();
    for (int i = 0; i < children.length; i++) {

        if (children[i].contains(str)) {
            System.out.println(children[i]);
            File fileremv = new File(pathout+"/"+children[i]);
            fileremv.delete();
        }
    }
    return true;
   }
}

  
