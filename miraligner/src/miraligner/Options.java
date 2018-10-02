/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package miraligner;
import com.beust.jcommander.Parameter;
import java.util.ArrayList;
import java.util.List;
/**
 *
 * @author lpantano
 */
public class Options {
  @Parameter
  public List<String> parameters = new ArrayList<String>();
 
  @Parameter(names = "-help", description = "help")
  public boolean help = false;
  
  @Parameter(names = "-v", description = "version")
  public boolean version = false;
 
  @Parameter(names = "-trim", description = "nt allowed for trimming processes")
  public String trim = "3";
  
  @Parameter(names = "-add", description = "nt allowed for addition processes")
  public String add = "3";
  
  @Parameter(names = "-sub", description = "nt allowed for substitution processes")
  public String sub = "1";
  
  @Parameter(names = "-s", description = "three letter code for species")
  public String species;
  
  @Parameter(names = "-db", description = "database")
  public String db;
  
  @Parameter(names = "-i", description = "input")
  public String input;
  
  @Parameter(names = "-o", description = "output")
  public String output;
  
  @Parameter(names = "-freq", description = "add freq information")
  public boolean freq = true;
  
  @Parameter(names = "-pre", description = "add reads mapping to precursor")
  public boolean pre;
  
  @Parameter(names = "-format", description = "format input")
  public String format = "none";

  @Parameter(names = "-minl", description = "minimum size")
  public Integer minl = 16;
  
}
