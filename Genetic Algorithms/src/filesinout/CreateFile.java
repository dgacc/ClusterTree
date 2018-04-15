 package filesinout;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.*;

import structures.Cluster;

public class CreateFile {
	
	public double[][] weightMatrix;
	double cost;
	private  ArrayList<Cluster> clusters= ReadFiles.clusters;
	private int num_vertex = ReadFiles.num_vertex;
	
	
	public  void writeFile(){
		
	 File newfile = new File("8i600NON.txt");
//	if(newfile.exists()){
//		System.out.println("this file already exists");
//	}else{
		try {
			newfile.createNewFile();
		}catch(Exception e){
			e.printStackTrace();
		}
		try {
			FileWriter instance = new FileWriter(newfile);
			BufferedWriter bufferedWriter = new BufferedWriter(instance);
			
			bufferedWriter.write("OPTIMAL COST : "+ cost +"\n");
			bufferedWriter.write("DIMENSION : "+ ReadFiles.num_vertex +"\n");
			bufferedWriter.write("NUMBER_OF_CLUSTERS: " + clusters.size()+"\n");
			
			for( int i = 0; i < num_vertex; i++){
//				bufferedWriter.write(i + 1 +" ");
				for( int j = 0; j < num_vertex; j++ ){
					bufferedWriter.write(weightMatrix[i][j]+" ");
				}
				bufferedWriter.write("\n");
			}
			
			bufferedWriter.write("CLUSTER_SECTION: \n");
			bufferedWriter.write("SOURCE_VERTEX: " +ReadFiles.root +"\n");
		
			for(int i = 0; i < clusters.size(); i ++){
				int numberClusterVertex = clusters.get(i).getCluster().size();
				bufferedWriter.write(i + 1 +" ");
				for(int j = 0 ; j < numberClusterVertex; j++){
					bufferedWriter.write(clusters.get(i).getCluster().get(j) +" ");
				}
				bufferedWriter.write(" -1");
				bufferedWriter.write("\n");
				
			}
			bufferedWriter.write("EOF");
			bufferedWriter.close();
			System.out.println("file written! ");
		}catch(Exception e){
			e.printStackTrace();
//		}
		
	}
	}
	

}
