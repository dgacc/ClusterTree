package clustertree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import filesinout.ReadFiles;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import operator.Selection;
import pruferpresentation.MFEAPruferNumber;
import random.MyRandom;
import structures.Individual;
import structures.Population;

public class CoevolutionCLT {
	private static double m_rate = 0.15;
	private static final double p_rate = 0.8; // crossover rate
	private static  int generation = 500;
	public static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();
	private static Evaluation evaluation = new Evaluation();
	
	
	
	public static void main(String arg[]){
		int besttour = 10 ;
		double[] popFitness;
		CoevolutionCLT coev = new CoevolutionCLT();
		MyRandom.setSeed(0);
		
		//initialize parameter
		int numberOfCluster = ReadFiles.clusters.size();
		int temp = (int)numberOfCluster/2;
		int [] startpoint = evaluation.getStartPoint(ReadFiles.clusters);
        int pointDivide	= startpoint[temp - 1];
        
        
        
        
		//initialize  the population 
		Population population = new Population();
		InitializeChromosome initializeChromosome = new InitializeChromosome();	
		population.population = initializeChromosome.initiaizePopulation(ReadFiles.clusters,populationLength, MyRandom.r);
		ArrayList<Population> listPopulation = new ArrayList<>();
		//evaluate initial population 
		popFitness = population.populationFitness1();
		// get the best tour
		int bestIndividual = population.getBestFitness(popFitness, populationLength);
		// divide into sub_populations 
		listPopulation = coev.dividePopulation(population.getIndividual(bestIndividual).getGene1(),
				evaluation.getStartPoint(ReadFiles.clusters), populationLength);
		// set active and frozen sub_population
	    Population  active = new Population();
	    Population  frozen = new Population(); 
	    active = listPopulation.get(0);
	    frozen = listPopulation.get(1);
	    active.setState(1);
	    active.setState(0);
		//combine  all of individual of active to the best of frozen
	    population = coev.combinePopulation(active, frozen.getIndividual(0), populationLength, pointDivide);
	    //evauuate Combine Pop;
	    popFitness = population.populationFitness1();
        //
		for(int i = 0; i< generation; i++){
			
			for(int j = 0; j < populationLength; j++ ){
				
				// select tour to do crossover to active sub_population
				Individual  father = new Individual();
				Individual  mother = new Individual();
				Individual  child1 = new Individual();
				Individual  child2 = new Individual();
				Individual  child3 = new Individual();
			    father = active.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));
			    mother = active.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));
			    child1.setGene(crossover.pruferNumberCrossover(father.getGene1(), mother.getGene1(),20,MyRandom.r).get(0));
			    child2.setGene(crossover.pruferNumberCrossover(father.getGene1(), mother.getGene1(),20,MyRandom.r).get(1));
			    //  apply mutation for frozen population 
			    child3.setGene(mutation.mutationPrufer(ReadFiles.clusters,	MyRandom.r,frozen.getIndividual(0).getGene1(),  startpoint, 0.03));
			}
			 population = coev.combinePopulation(active, frozen.getIndividual(0), populationLength, pointDivide);
			 popFitness = population.populationFitness1();
			 bestIndividual = population.getBestFitness(popFitness, populationLength);
			 
		}
	}
	
	
	/**
	 * divide population into 2 subPopulation
	 */
	public ArrayList<Population> dividePopulation(int[] pruferCode, int[] startpont, int subPopLength){
	    ArrayList<Population> twoPopulation = new  ArrayList<Population>();
	    InitializeChromosome initializeChromosome = new InitializeChromosome();
		int numberOfClusters = ReadFiles.clusters.size();
		Population subPopulation1 = null;
		Population subPopulation2 = null;
		Individual indv1 = new Individual();
		Individual indv2 = new Individual();
		int[] pruferCode1 = null;
		int[] pruferCode2 = null;
		int geneLength =  0; 
		int sub_NumberClt = (int)numberOfClusters/2;
		int pointDevide = startpont[sub_NumberClt - 1];
		
		for( int i = 0 ; i < subPopLength; i++){
			
			ArrayList<Integer> temp  = new ArrayList<>();
			for(int j = 0; j < numberOfClusters - 1; j++){
				
				if(j == 0){
					ArrayList<Integer> temp1 = new ArrayList<>();
					for(int k = 0; k < startpont[0]; k++){
						temp1.add(pruferCode[k]);
						}
					Collections.shuffle(temp1);
		          temp.addAll(temp1);
				}else{
					ArrayList<Integer> temp1 = new ArrayList<>();
				for(int k = startpont[j -1]; k < startpont[j]; k++){
					temp1.add(pruferCode[k + startpont[j]]);
				}
				Collections.shuffle(temp1);
				temp.addAll(temp1);
				}
				}
			
			for(int j = 0; j< geneLength; j++){
				if( j <  pointDevide){
					pruferCode1[j] = temp.get(j);
				}else{
					if(j <  geneLength - numberOfClusters){
					pruferCode2[j - pointDevide] = temp.get(j);
				    }else{
					pruferCode2[j] = pruferCode[j];
				}
			  }
			}	
			
			indv1.setGene(pruferCode1);
			indv2.setGene(pruferCode2);
			subPopulation1.addIndiv(indv1);
			subPopulation2.addIndiv(indv2);
		}
		twoPopulation.add(subPopulation1);
		twoPopulation.add(subPopulation2);
		return twoPopulation;
	}
	
	/**
	 * set active and frozen to the sub-population 
	 */
	public void setStatePopulation(){
     
		
	}
	/**
	 * match randomly two individual in the first generation
	 */
	public Population matchSubPopulation(Population subPop1, Population subPop2, int popSize, Random r){
		MFEAPruferNumber mfeaPruferNumber  = new MFEAPruferNumber();
		int genleng1 = mfeaPruferNumber.calCulateGenLength(ReadFiles.clusters);
		Population population  = null;
		Individual individual;
		for(int i = 0; i < genleng1; i++){
			for(int  j = 0 ; j< genleng1; j++){
				
			}	
		}
   return population;
		
	}
	/**
	 * combine the best individual in frozen sub_population with all of individual in Active sub-population
	 */
	public Population combinePopulation( Population active , Individual frozen, int popLength, int pointDivide){
    
		Population population = new Population();
		Individual individual = new Individual();
//		double[] frozenFitness = frozen.populationFitness1();
//		int bestIndividual  = 0;
		int geneLength  = 0 ;// chi dung tap thoi de debug 
//		
//		// get the best individual in frozen population
//		double bestFzFiness = frozenFitness[0]; 
//		for(int i = 0; i < popLength; i++){
//	         if(bestFzFiness < frozenFitness[i]){
//	        	 bestIndividual = i;
//	        	 bestFzFiness = frozenFitness[i];
//	         }
//		}
		int[] prozenChro = frozen.getGene1();
		// match the best individual in frozen Sub_Population to all of individual in active  Sub_population
		for(int i = 0; i < popLength; i++){
			int[] pruferCode = null;
			int[] activeChro = active.getIndividual(i).getGene1();
		
			for(int j = 0 ; j < geneLength ; j++){
				if( j < pointDivide){
					pruferCode[j] =   activeChro[j];
				}else{
					pruferCode[j] =  prozenChro[j - pointDivide];
				}
			}
			individual.setGene(pruferCode);
			population.addIndiv(individual);
     	}
		return population;
		
	}
	/**
	 * calculate the population, then calculate the best fitness
	 * @param population
	 * @return the index of fitness
	 */
	public int  bestTour(Population population){
		double[] popFitness = population.populationFitness1();
		//get best fitness
		double bestFitness  = popFitness[0];
		int bestFitnessIndex = 0;
		for(int i = 0 ; i< populationLength; i++){
			if(popFitness[i] < bestFitness){
				bestFitness = popFitness[i];
				bestFitnessIndex = i;
			}
		}
		return bestFitnessIndex;
	}
	
	

}
