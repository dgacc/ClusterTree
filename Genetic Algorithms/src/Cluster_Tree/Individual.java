package Cluster_Tree;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Individual {
	public  static int num_vertex = ReadFiles.num_vertex;
	public  double[][] gene =  new double[num_vertex][num_vertex]; 
	InitializeChromosome initializechromosome = new InitializeChromosome();
	Evaluation evaluation = new Evaluation();
	
	
	public Individual() {
	}
    public void inintilizeIndividual(){
    	gene = initializechromosome.clusterPrimRST(ReadFiles.weightMatrix, ReadFiles.clusters, num_vertex);
 
    }
    // getter
    public double[][] getGene(){
    	return gene;
    }
    //setter
    public void setGene(double[][] gene){
    	 this.gene = gene;	
    }
    
    // get fitness
    public double getFitness(Individual individual){
    	 return evaluation.evaluation(ReadFiles.weightMatrix, individual.getGene(), num_vertex, ReadFiles.root);
    }
   
}