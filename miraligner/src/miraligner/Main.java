/*

 */

package miraligner;

import java.io.FileNotFoundException;
import java.io.IOException;
import com.beust.jcommander.*;

public class Main {


    public static void main(String[] args) throws FileNotFoundException, IOException {
        String format = "None";
        String test="notest";
        if (test.equals("test")){
            map.readseq("test/test.fa","DB","hsa",1,3,3,"fasta","test/test",false,true,16);
            System.exit(0);
        }
        Options jct = new Options();
        JCommander jc = new JCommander(jct, args);
        if (jct.help | args.length<4 ){
            jc.usage();
            System.out.println("\njava -jar miraligner.jar -minl 16 -sub mismatches -trim trimming-nts -add addition-nts -s species -i read_seq_file -db miRBase_folder_files -o output_file");
            System.out.println("\nexample:java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s hsa -i test/test.fa -db DB -o test/out");
            System.out.println("example: see output at miraligner/test/output.mirna & miraligner/test/output.mirna.opt");
            System.out.println("\n");
            System.exit(0);
        }
        if (jct.version) {
            System.out.println("version 2");
            System.exit(0);
        }
        //check input file
        if ("none".equals(jct.format)){
            boolean f=tools.checkinput(jct.input);
            boolean ftab=tools.checkinputtab(jct.input);
            if (f){
                format="fasta";
            }else if (ftab){
                format="tab";
            }else{
                System.err.println("no format file recognized (fasta or tabular)");
                System.exit(1);
            }
        }else{
            format = jct.format;
        }
        //check directory
        
        //check species
        boolean sp=tools.checksp(jct.db,jct.species);
        int mism=Integer.parseInt(jct.sub);
        if (mism>1){
           System.out.println("Only allowed 0/1 mismatch");
           sp=false;
        }
        int trim=Integer.parseInt(jct.trim);
        if (trim>3){
           System.out.println("Only allowed <=3 nucleotides as trimming");
           sp=false;
        }
        int add=Integer.parseInt(jct.add);
        if (add>3){
           System.out.println("Only allowed <=3 nucleotides as addition");
           sp=false;
        }
        if (jct.minl<16){
           System.out.println("Only allowed >=16 minimum size");
           jct.minl=16;
        }
        if (sp  ){
            System.out.println("Go to mapping...");
            System.out.println("Mismatches: "+jct.sub);
            System.out.println("Trimming: "+jct.trim);
            System.out.println("Addition: "+jct.add);
            System.out.println("Species: "+jct.species);
            
            map.readseq(jct.input,jct.db,jct.species,mism,trim,add,format,jct.output,jct.freq,jct.pre,jct.minl);

       }else{
           System.err.println("species not found: "+jct.species);
           System.exit(1); 
        }
    }



}
