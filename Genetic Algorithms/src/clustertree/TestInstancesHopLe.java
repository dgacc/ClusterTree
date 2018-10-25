package clustertree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import filesinout.ReadFiles;
import operator.InitializeChromosome;

public class TestInstancesHopLe {

	public static void main(String[] args) {

		String[] test = {"100i1000-410","100i1500-506","100i2000-604","100i2500-708","100i3000-803","10i1000-407","10i1500-503","10i300-109","10i400-206","10i500-305","150i1000-411","150i1500-507","150i2000-605","150i2500-709","150i3000-804","15i300-110","15i400-207","15i500-306","200i2000-606","200i2500-710","200i3000-805","20i1000-408","20i1500-504","20i2000-602","20i2500-706","20i300-111","20i3000-801","20i400-208","20i500-307","25i300-112","25i400-209","25i500-308","50i1000-409","50i1500-505","50i2000-603","50i2500-707","50i3000-802","5i300-108","5i400-205","5i500-304"};

		 String dirname = "InstanceHopLe";
	      File d = new File(dirname);
	      // Bay gio tao thu muc.
	      d.mkdirs();
	      
	     PrintWriter pw = null;
	     try {

				pw = new PrintWriter(new FileWriter(dirname+"/"+"Type5Large.hl", false));

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

       		pw.print( test[k]+" ");

       		count++;
       	}
	}pw.close();
	System.out.println("Done!");
	System.out.println("Number Of Instances : " +test.length);
	System.out.println("Number Of Instances  hop le : " + count);
	
	
}
	
	}
