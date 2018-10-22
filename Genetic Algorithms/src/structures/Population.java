package structures;

import java.util.ArrayList;

import filesinout.ReadFiles;
import operator.Evaluation;

public class Population {
	public static int populationLength = 100; // default length of population
	public ArrayList<Individual> population = new ArrayList<Individual>();
    Individual individual = new Individual();
    Evaluation evaluation = new Evaluation();
    private int state; // 1 = activate and 0 = frozen
    
	
	public Population() {
	};
	
	// get gene of an Individual
	public Individual getIndividual(int i) {
		return population.get(i);
	}
	
	public void addIndiv(Individual individual){
		population.add(individual);
	}
	
	public ArrayList<Individual> getPopulation(){
		return this.population;
	}

	
	public double[] populationFitness(){
		double[] popFitness = new double[populationLength];
		for( int i = 0; i < populationLength; i++){
			double[][] tree = population.get(i).getGene();
			popFitness[i] = evaluation.evaluation(ReadFiles.weightMatrix, tree, ReadFiles.num_vertex, ReadFiles.root);
		}
		return popFitness;
	}
	/**
	 * evaluate the fitness for population, which use prufer number presentation
	 * @return the array of fitness of  the population 
	 */
			
	public double[] populationFitness1(){
		double[] popFitness = new double[populationLength];
		for( int i = 0; i < populationLength; i++){
			double[][] tree = evaluation.decodingPruferNumber(population.get(i).getGene1(), evaluation.getStartPoint(ReadFiles.clusters),
					ReadFiles.clusters, ReadFiles.num_vertex);
			
			popFitness[i] = evaluation.evaluation(ReadFiles.weightMatrix, tree, ReadFiles.num_vertex, ReadFiles.root);
		}
		return popFitness;
	}
	
	/***
	 * Get index of the best fitness of  individual in Population  
	 * @param popFitness : the double[] fitness of each individual in population 
	 * @param popLength :  the length of population
	 * @return the index of the best fitness
	 */
	public int getBestFitness(double[] popFitness, int popLength){
		int bestIndividual = 0;
		double bestFiness = popFitness[0];
		for(int i = 0 ; i < popLength; i++){
			if(popFitness[i] < bestFiness){
				bestIndividual = 0;
				bestFiness = popFitness[0];
			}
		}
		return bestIndividual;	
	}
	
	
	public int getState(){
	    return this.state;
	}
	public void setState(int state){
		this.state = state;
	}
	
}
	
