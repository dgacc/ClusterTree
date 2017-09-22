package Cluster_Tree;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Individual {
	public int num_vertex = ReadFiles.num_vertex;
	public double[][] gene =  new double[num_vertex][num_vertex]; 
	InitializeChromosome chromosome = new InitializeChromosome();
	
	public Individual() {
	}
    public void inintilizeIndividual(){
    	gene = chromosome.clusterPrimRST(ReadFiles.distances, ReadFiles.clusters, num_vertex);
    	
    }
    // getter
    public double[][] getGene(){
    	return gene;
    }
    
   
	/*public void printIndiv() {
		for (int i = 0; i < num_vertex -2; i++) {
			System.out.print(gene.get(i) + " ");
		}
		System.out.println();
		System.out.print("The Spanning  Tree  Cost is :  ");
		System.out.print(this.getFitness());
		System.out.println();
		System.out.println("________________________________________________________________________________________________________________________________________________");

		System.out.println();
	}
	*/
}