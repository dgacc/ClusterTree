package clustertree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import filesinout.ReadFiles;
import operator.InitializeChromosome;

public class TestInstancesHopLe {

	public static void main(String[] args) {
		String[] test = {"6i300","6i350","6i400","6i450"};
		 String dirname = "InstanceHopLe";
	      File d = new File(dirname);
	      // Bay gio tao thu muc.
	      d.mkdirs();
	      
	     PrintWriter pw = null;
	     try {
				pw = new PrintWriter(new FileWriter(dirname+"/"+"Type3_large.hl", false));
			} catch (IOException e) {
				e.printStackTrace();
			}
	     int count = 0;
		for(int k= 0; k < test.length; k++){
			 boolean flag = true;
					ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/"+ test[k]+".clt");
	      if(ReadFiles.clusters.size() < 3){
	    	  flag = false;
	      }
       	for( int i = 0; i < ReadFiles.clusters.size(); i++){
       		if(ReadFiles.clusters.get(i).getCluster().size() < 3){
       			flag = false;
       			break;
       		}    		
       	}
       	if(flag){
       		pw.print("\""+ test[k]+"\",");
       		count++;
       	}
	}pw.close();
	System.out.println("Done!");
	System.out.println("Number Of Instances : " +test.length);
	System.out.println("Number Of Instances  hop le : " + count);
	
	
}
	
	}
