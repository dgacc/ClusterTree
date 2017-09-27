package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;

public class Population {
	public static int populationLength = 100; // default length of population
	public ArrayList<Individual> population = new ArrayList<Individual>(populationLength);
    public Evaluation evaluation = new Evaluation();
    public InitializeChromosome initializeChromosome = new InitializeChromosome();
    
    
	public Population() {
	};
	
	// get gene of an Individual
	public Individual getIndividual(int i) {
		return population.get(i);
	}
	
	public void addIndiv(Individual individual){
		population.add(individual);
	}
	public  double[][] buildWeightMatrixIndividual(Individual individual){
		double[][] weightMatrix = ReadFiles.weightMatrix ;
		double[][] weightMatrix1 = individual.getGene();
		double[][] weightMatrix2 = new double[ReadFiles.num_vertex][ReadFiles.num_vertex];
		int num_vertex = ReadFiles.num_vertex;
		for(int i = 0; i < num_vertex; i ++){
			 for( int j = 0; j < num_vertex; j ++ ){
				 if(weightMatrix1[i][j] > 0)
				 weightMatrix2[i][j] = weightMatrix[i][j];
			 }
		 }
		return weightMatrix2;
	}
	/*public double[][] getBestIndividual(){
		double bestFitness = 0;
		for(int i = 0 ; i < populationLength; i ++){
			evaluation.clusterEvaluate(weightMatrix, edgeMatrix, num_vertex, startVertex)
		}
		
	}*/

	
	}
