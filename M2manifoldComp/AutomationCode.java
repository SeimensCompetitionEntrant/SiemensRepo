import java.util.Scanner;
import java.util.ArrayList;
import java.io.*;
import java.text.*;

public class AutomationCode{
  private String necPath;
  private String folderPath;
  private String simulFileName;
  private int numAnts;
  private int pattern;
  private double freq;
  private double wvl;
  private double spacing;
  private double dist;
  private final double c;


  public static void main(String[] args) {
    AutomationCode ac =  new AutomationCode();
    ac.requestAndCalc();
    ac.generateInputFiles();
    ac.generateOutputFilesandMove();
    ac.callOctave();
  }

  public AutomationCode(){
    necPath = new String("");
    File tempFile = new File(".");
    folderPath = tempFile.getAbsolutePath();
    folderPath = folderPath.substring(0,folderPath.length()-2);
    simulFileName = "";
    numAnts = 0;
    pattern = 0;
    freq = 0.0;
    wvl = 0.0;
    spacing = 0.0;
    dist = 0.0;
    c = 299792458.0;
  }

  public void requestAndCalc(){
    Scanner input = new Scanner(System.in);
    necPath = Prompt.getString("Please enter the full path to the 4nec2 installation. Please make sure that the installation of 4nec2 and this folder are on the same drive: ");
    numAnts = Prompt.getInt("\nNumber of antennas: ");
    pattern = Prompt.getInt("\nPattern? 1 for ULA, 2 for UCA: ");
    freq = Prompt.getDouble("\nCarrier frequency, in MHz: ");
    spacing = Prompt.getDouble("\nSpacing between antennas in wavelengths: ");
    wvl = c/(freq * 1000000.0);
    dist = spacing*wvl;
  }

  public void generateInputFiles(){
    DecimalFormat df = new DecimalFormat("00.00000");
    File f = null;
    PrintWriter outFile = null;
    String antFileName = "";
    for(int i = 0; i < numAnts; i++){
      antFileName = "ant" + (i+1);
      simulFileName += antFileName;
      f = new File("./workingCouplingComp/NECuserMadeFiles/vDipoleArray/test/" + antFileName + ".nec");
      try{
  			outFile = new PrintWriter(f);
  		}catch(Exception e){
  			System.out.println("\nSorry, but the file could not be created.\n");
  			System.exit(2);
  		}
      outFile.print("CM\r\nCE\r\n");
      for(int k = 0; k < numAnts; k++){
        outFile.print("GW  " + (k+1) + "  9 0 " + df.format(k*dist) + " " + df.format(-wvl/4) + " 0 " + df.format(k*dist) + " " + df.format(wvl/4) + " .0001\r\n");
      }
      outFile.print("GE  0\r\n");
      for(int k = 0; k < numAnts; k++){
        outFile.print("LD  4 " + (k+1) + " 5 5 50\r\n");
      }
      outFile.print("GN -1\r\n");
      outFile.print("EK\r\n");
      outFile.print("EX  0 " + (i+1) + " 5 0 1 0 0\r\n");
      outFile.print("FR  0 0 0 0 " + freq + " 0\r\n");
      outFile.print("RP  0 37  73    1003 -180     0         5         5\r\n");
      outFile.print("EN");
      outFile.close();
    }
    f = new File("./workingCouplingComp/NECuserMadeFiles/vDipoleArray/test/" + simulFileName + ".nec");
    try{
      outFile = new PrintWriter(f);
    }catch(Exception e){
      System.out.println("\nSorry, but the file could not be created.\n");
      System.exit(2);
    }
    outFile.print("CM\r\nCE\r\n");
    for(int k = 0; k < numAnts; k++){
      outFile.print("GW  " + (k+1) + "  9 0 " + df.format(k*dist) + " " + df.format(-wvl/4) + " 0 " + df.format(k*dist) + " " + df.format(wvl/4) + " .0001\r\n");
    }
    outFile.print("GE  0\r\n");
    for(int k = 0; k < numAnts; k++){
      outFile.print("LD  4 " + (k+1) + " 5 5 50\r\n");
    }
    outFile.print("GN -1\r\n");
    outFile.print("EK\r\n");
    for(int k = 0; k < numAnts; k++){
      outFile.print("EX  0 " + (k+1) + " 5 0 1 0 0\r\n");
    }
    outFile.print("FR  0 0 0 0 " + freq + " 0\r\n");
    outFile.print("RP  0 37  73    1003 -180     0         5         5\r\n");
    outFile.print("EN");
    outFile.close();
  }

  public void generateOutputFilesandMove(){
    for(int i = 0; i < numAnts; i++){
      String antFileName = "ant" + (i+1);
      try {
        String s = null;
        Process p = Runtime.getRuntime().exec("4nec2.exe " + folderPath + "\\workingCouplingComp\\NECuserMadeFiles\\vDipoleArray\\test\\" + antFileName + ".nec -I");

        BufferedReader stdInput = new BufferedReader(new
            InputStreamReader(p.getInputStream()));

        BufferedReader stdError = new BufferedReader(new
            InputStreamReader(p.getErrorStream()));

        // read the output from the command
        System.out.println("Here is the standard output of the command:\n");
        while ((s = stdInput.readLine()) != null) {
          System.out.println(s);
        }

        // read any errors from the attempted command
        System.out.println("Here is the standard error of the command (if any):\n");
        while ((s = stdError.readLine()) != null) {
          System.out.println(s);
        }
      }
      catch (IOException e) {
        System.out.println("exception happened - here's what I know: ");
        e.printStackTrace();
        System.exit(-1);
      }
      File f = new File("./workingCouplingComp/NECuserMadeFiles/vDipoleArray/test/" + antFileName + ".txt");
      if(f.exists()){
        f.delete();
      }
      f = new File(necPath + "/out/" + antFileName + ".out");
      f.renameTo(new File("./workingCouplingComp/NECuserMadeFiles/vDipoleArray/test/" + antFileName + ".txt"));
    }
    try{
      String s = null;
      Process p = Runtime.getRuntime().exec("4nec2.exe " + folderPath + "\\workingCouplingComp\\NECuserMadeFiles\\vDipoleArray\\test\\" + simulFileName + ".nec -I");

      BufferedReader stdInput = new BufferedReader(new
          InputStreamReader(p.getInputStream()));

      BufferedReader stdError = new BufferedReader(new
          InputStreamReader(p.getErrorStream()));

      // read the output from the command
      System.out.println("Here is the standard output of the command:\n");
      while ((s = stdInput.readLine()) != null) {
        System.out.println(s);
      }

      // read any errors from the attempted command
      System.out.println("Here is the standard error of the command (if any):\n");
      while ((s = stdError.readLine()) != null) {
      System.out.println(s);
      }
    }
    catch (IOException e) {
      System.out.println("exception happened - here's what I know: ");
      e.printStackTrace();
      System.exit(-1);
    }
    File f = new File("./workingCouplingComp/NECuserMadeFiles/vDipoleArray/test/" + simulFileName + ".txt");
    if(f.exists()){
      f.delete();
    }
    f = new File(necPath + "/out/" + simulFileName + ".out");
    f.renameTo(new File("./workingCouplingComp/NECuserMadeFiles/vDipoleArray/test/" + simulFileName + ".txt"));
  }


  public void callOctave(){
    try {
      String s = null;
      System.out.println("octave " + folderPath + "\\workingCouplingComp\\MatlabFiles\\M2vFinalArray.m");
      Process p = Runtime.getRuntime().exec("octave " + folderPath + "\\workingCouplingComp\\MatlabFiles\\M2vArrayFinal.m");

      BufferedReader stdInput = new BufferedReader(new
          InputStreamReader(p.getInputStream()));

      BufferedReader stdError = new BufferedReader(new
          InputStreamReader(p.getErrorStream()));

      // read the output from the command
      System.out.println("Here is the standard output of the command:\n");
      while ((s = stdInput.readLine()) != null) {
        System.out.println(s);
      }

      // read any errors from the attempted command
      System.out.println("Here is the standard error of the command (if any):\n");
      while ((s = stdError.readLine()) != null) {
        System.out.println(s);
      }
    }
    catch (IOException e) {
      System.out.println("exception happened - here's what I know: ");
      e.printStackTrace();
      System.exit(-1);
    }
  }
}
