package pruferpresentation;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import dislay.Paint;
import dislay.Windows;
import edgepresentation.MEF2Instances;
import filesinout.ReadFiles;
import operator.ChromosomeCmp;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import random.MyRandom;
//"10i30-17","10i30-17","10i45-18","10i45-18","10i60-21","10i60-21","10i65-21","10i65-21","10i70-21","10i70-21","10i75-22","10i75-22","10i90-33","10i90-33"};
//"4i200a","4i200h","4i200x1","4i200x2","4i200z","4i400a","4i400h","4i400x1","4i400x2","4i400z"
//"10berlin52","10eil51","10eil76","10kroB100","10pr76","10rat99","10st70","15berlin52" type 1
//"10berlin52","10eil51","10eil76","10kroB100","10pr76","10rat99",
//"15eil51","15eil76","15pr76","15st70","25eil101","25kroA100","25lin105","25rat99","50eil101","50kroA100",
//"50kroB100","50lin105","50rat99","5berlin52","5eil51","5eil76","5pr76","5st70"
import structures.Cluster;
import structures.Individual;


public class MFEAPruferNumber {
	ReadFiles re = new ReadFiles();
	private Evaluation eva = new Evaluation();
	private Individual chromo = new Individual();
	private  Mutations mutation = new Mutations();
	private Crossover crossover = new Crossover();
	private static double rmp = 0.5;
	public int maxGroup;  
	public static int defaultPopLength = 100;

	public static void main(String[] args) {
		MFEAPruferNumber mfeaPruferNumber = new MFEAPruferNumber();
		String[] test = {"5i500-304","6i500","7i70-21","9pr439-3x3","2lin105-2x1","7i70-21","6i300","6i350","6i400","6i450"};
		for(int k= 0; k < test.length; k= k+2){
//		ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/"+ test[k]+".clt");
//		ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/"+ test[k]+".clt");
//		ReadFiles.clusterReadFiles1("C:/Users/TrungTB/Desktop/test/"+test[k+1]+".clt");
		ReadFiles.clusterReadFiles("test/"+ test[k]+".clt");
		ReadFiles.clusterReadFiles1("test/"+test[k+1]+".clt");
		FileOutputStream thu = null;
		System.out.println("-------------------------------------------------------------------------------");		
		// write file to use later
		 String dirname = "Results/"+test[k]+"_and_"+test[k+1];
	      File d = new File(dirname);
	      // Bay gio tao thu muc.
	      d.mkdirs();
		for( int loop = 0; loop < 30; loop++){
			 PrintWriter pw = null;
			 PrintWriter pw1 = null;
			 PrintWriter pw3 = null;
//			 PrintWriter pw4 = null;
			
			
				try {
					pw = new PrintWriter(new FileWriter(dirname+"/"+test[k]+"-seed("+loop+").opt", true));
					pw1 = new PrintWriter(new FileWriter(dirname+"/"+test[k+1]+"-seed("+loop+").opt", true));
					pw3 = new PrintWriter(new FileWriter(dirname+"/"+test[k]+"_and_"+test[k+1]+"_seed("+loop+").gen", true));
				} catch (IOException e) {
					e.printStackTrace();
				}
				pw3.println("Generations \t"+test[k]+"\t\t" +test[k+1]);
        long start = System.currentTimeMillis();
		double[] bestCost =  new double[2];
		ArrayList<Individual> bestIndividual = new ArrayList<Individual>(); 
		MyRandom.setSeed(loop);
		MEF2Instances mef2 = new MEF2Instances();
		InitializeChromosome initChromo = new InitializeChromosome();
		initChromo.clutersVerticesInf(ReadFiles.clusters, ReadFiles.clusters1);
		initChromo.buildCluster();
//		double muRate = mfeaPruferNumber.calMuRate(InitializeChromosome.maxClusters);
		double muRate  = 0.05;
		ArrayList<Individual> pop = new ArrayList<Individual>();
		Individual individual = new Individual();
		// ininitlize new population
	    pop  = mfeaPruferNumber.initiaizePopulation(InitializeChromosome.maxClusters,defaultPopLength, MyRandom.r);
		
		calculate(pop,defaultPopLength);
		calculateFactorialRank(pop, defaultPopLength);
		calculateScalarFitness(pop, defaultPopLength);
		
	
		sortByCostIndex(pop,defaultPopLength, 0);
		bestCost[0] = pop.get(0).getFactorialCost()[0];
//		System.out.print(pop.get(0).getFactorialCost()[0] + " ");
		bestIndividual.add(0, pop.get(0));
		sortByCostIndex(pop,defaultPopLength, 1);
		bestCost[1] = pop.get(0).getFactorialCost()[1];
//		System.out.println(pop.get(0).getFactorialCost()[1]);
		bestIndividual.add(1, pop.get(0));
		
//		Windows windows1 = new Windows();
//		Windows windows = new Windows();
//		windows.runWindow(" cac cay");
		for( int i = 0; i < 500; i++){
			 ArrayList<Individual> newPop = new ArrayList<Individual>();
			ArrayList<Individual> tempPop = new ArrayList<Individual>();
			  while(tempPop.size() < defaultPopLength){
				    double r = MyRandom.r.nextDouble();
					int par1 = MyRandom.r.nextInt(defaultPopLength);
					int par2 = MyRandom.r.nextInt(defaultPopLength);
					if (pop.get(par1).getSkillFactor() == pop.get(par2).getSkillFactor() || r < rmp)
					{
						
						Individual child = new Individual();
						Individual child1 = new Individual();
						child.setGene(mfeaPruferNumber.crossover.pruferNumberCrossover(pop.get(par1).getGene1(), pop.get(par1).getGene1(),
								mfeaPruferNumber.calCulateGenLength(InitializeChromosome.maxClusters), MyRandom.r).get(0));
						child1.setGene(mfeaPruferNumber.crossover.pruferNumberCrossover(pop.get(par1).getGene1(), pop.get(par1).getGene1(),
								mfeaPruferNumber.calCulateGenLength(InitializeChromosome.maxClusters) , MyRandom.r).get(1));
//						double rfs1 =  0 + (1 - 0) * MyRandom.r.nextDouble();
						if (MyRandom.r.nextDouble() < 0.5)
                        {
						 child.setSkillFactor( pop.get(par1).getSkillFactor());
                        }
                        else
                        {
                        	child.setSkillFactor( pop.get(par2).getSkillFactor());
                        }
						if (MyRandom.r.nextDouble() < 0.5)
                        {
							child1.setSkillFactor( pop.get(par1).getSkillFactor());
                        }
                        else
                        {
                        	child1.setSkillFactor( pop.get(par2).getSkillFactor());
                        }
						
						tempPop.add(child);
						tempPop.add(child1);
						
					
					}else{
						
						pop.get(par1).setGene(mfeaPruferNumber.mutation.mutationPrufer(InitializeChromosome.maxClusters,
								MyRandom.r, pop.get(par1).getGene1(), mfeaPruferNumber.eva.getStartPoint(InitializeChromosome.maxClusters), muRate));
						tempPop.add(pop.get(par1));
						pop.get(par1).setGene(mfeaPruferNumber.mutation.mutationPrufer(InitializeChromosome.maxClusters,
								MyRandom.r, pop.get(par2).getGene1(), mfeaPruferNumber.eva.getStartPoint(InitializeChromosome.maxClusters), muRate));
						tempPop.add(pop.get(par2));
					}  
			  }
			  if(tempPop.size() > defaultPopLength){
				  tempPop.remove(defaultPopLength);
			  }
			  calculate(tempPop, defaultPopLength);
			  
			
				  newPop.addAll(tempPop);
				  newPop.addAll(pop);
				  
			  
			  calculateFactorialRank(newPop, defaultPopLength*2);
			  
			  
			  
			  sortByCostIndex(newPop,defaultPopLength*2, 0);
			  if(newPop.get(0).getFactorialCost()[0] < bestCost[0]){
				  bestCost[0] = newPop.get(0).getFactorialCost()[0];
				  bestIndividual.set(0, newPop.get(0));
			  }
			  pw3.print(i+"\t"+bestCost[0] + "\t");
//			  System.out.print(newPop.get(0).getFactorialCost()[0] + " ");
			  sortByCostIndex(newPop, defaultPopLength*2, 1);
			  if(newPop.get(0).getFactorialCost()[1] < bestCost[1]){
				  bestCost[1] = newPop.get(0).getFactorialCost()[1];
				  bestIndividual.set(1, newPop.get(0));
			  }
			  pw3.println(bestCost[1]);
			  
//			  Paint  p = new Paint();
//			  double[][] caytest = mfeaPruferNumber.eva.decodingPruferNumber(newPop.get(0).getGene1(), mfeaPruferNumber.eva.getStartPoint(InitializeChromosome.maxClusters),
//		        		ReadFiles.clusters1,InitializeChromosome.maxClusters, ReadFiles.num_vertex1 , MyRandom.r);
//			  
//			  p.setPaint(caytest, ReadFiles.vertices1, ReadFiles.clusters1, ReadFiles.num_vertex1, 
//					  tempPop.get(0).getFactorialCost()[1], 0, ReadFiles.root1);
//			  windows.addPaint(p);
//			 
			  
			  calculateScalarFitness(newPop, defaultPopLength*2);
			 
				for (int j = 0; j < defaultPopLength; j++) {
					pop.set(j,newPop.get(j));
					
				}
			
		}
		pw3.close();
		
		long end = System.currentTimeMillis();
		System.out.println("|"+loop+"\t|"+test[k]+"\t|"+test[k+1]+"\t|"+bestCost[0]+"\t|"+bestCost[1]+"\t|"+getDateFromMillis(end - start)+"|");
		 pw.println("Name: "+test[k]);
		 pw.println("Seed: " + loop);
		 pw.println("Fitness: " +bestCost[0]);
		 pw.println("Time: " + getDateFromMillis(end - start));
//		 for(int m = 0; m < ReadFiles.num_vertex; m ++){
//			 pw.println();
//			 for( int n = 0; n < ReadFiles.num_vertex; n++){
//				 if(  bestIndividual.get(0).getGene1()[m][n] > 0){
//					 bestIndividual.get(0).getGene1()[m][n] = 1.0f; }
//				 pw.print((int)bestIndividual.get(0).getGene1()[m][n]+"\t");
//				
//			 }
//		 }
		 pw.close();
//		 windows1.runWindow("Cay tot nhat");
//		 Paint  p = new Paint();
//		  double[][] caytest = mfeaPruferNumber.eva.decodingPruferNumber(bestIndividual.get(0).getGene1(), mfeaPruferNumber.eva.getStartPoint(InitializeChromosome.maxClusters),
//	        		ReadFiles.clusters1, InitializeChromosome.maxClusters, ReadFiles.num_vertex1 , MyRandom.r);
//		  
//		  p.setPaint(caytest, ReadFiles.vertices1, ReadFiles.clusters1, ReadFiles.num_vertex1, 
//				  bestIndividual.get(0).getFactorialCost()[1], 0, ReadFiles.root1);
//		  windows1.addPaint(p);
		 
		 pw1.println("Name: "+test[k+1]);
		 pw1.println("Seed: " + loop);
		 pw1.println("Fitness: " +bestCost[1]);
		 pw1.println("Time: " + getDateFromMillis(end - start));
//		 for(int m = 0; m < ReadFiles.num_vertex1; m ++){
//			 pw1.println();
//			 for( int n = 0; n < ReadFiles.num_vertex1; n++){
//				 if(  bestIndividual.get(0).getGene1()[m][n] > 0){
//					 bestIndividual.get(0).getGene1()[m][n] = 1.0f; }
//				 pw1.print((int)bestIndividual.get(0).getGene1()[m][n]+"\t");
//				
//			 }
//		 }
		 pw1.close();
		 
	  }
	 
	
		}

	
	}
	
	
	public static void calculateFactorialRank(ArrayList<Individual> pop, int pop_size) {
		sortByCostIndex(pop,pop_size, 0);
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).setFactorialRank(i + 1, 0);
//			 System.out.println(pop.get(i).getFactorialCost()[0] + " ");
		}
		sortByCostIndex(pop,pop_size, 1);
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).setFactorialRank(i + 1, 1);
		}
	}
	
	public static void calculateScalarFitness(ArrayList<Individual> pop, int pop_size){
		for (int i = 0; i < pop_size; i++) {
			if (pop.get(i).getFactorialRank()[0] < pop.get(i).getFactorialRank()[1]) {
				pop.get(i).setSkillFactor(0);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[0]) + 1);
			} else {
				pop.get(i).setSkillFactor(1);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[1]) + 1);
			}
		}
		Collections.sort(pop, ChromosomeCmp.compareByScalarFitness);
	}

	
	/**
	 * calculate  factorial cost
	 * @param pop
	 */
	public static void calculate(ArrayList<Individual> pop, int pop_size) {
		MEF2Instances mc = new MEF2Instances();
		InitializeChromosome initChromo = new InitializeChromosome();
		double weightMatrix[][] = ReadFiles.weightMatrix;
		double[][] weightMatrix1 = ReadFiles.weightMatrix1;
		int startVertex1 = ReadFiles.root;
		int startVertex2 = ReadFiles.root1;
		int num_vertex1 = ReadFiles.num_vertex;
		int num_vertex2 = ReadFiles.num_vertex1;
		
		
		
		for (int i = 0; i < pop_size; i++) {
			double[] temp = new double[2];
			  double[][] TreeForEachInst;
			  double[][] TreeForEachInst1;
		        TreeForEachInst = mc.eva.decodingPruferNumber(pop.get(i).getGene1(), mc.eva.getStartPoint(InitializeChromosome.maxClusters),
		        		ReadFiles.clusters, InitializeChromosome.maxClusters, num_vertex1 , MyRandom.r);
		        TreeForEachInst1 = mc.eva.decodingPruferNumber(pop.get(i).getGene1(), mc.eva.getStartPoint(InitializeChromosome.maxClusters),
		        		ReadFiles.clusters1,  InitializeChromosome.maxClusters, num_vertex2 , MyRandom.r);
//			
			temp[0] = mc.eva.evaluation( weightMatrix, TreeForEachInst, num_vertex1, startVertex1);
			temp[1] = mc.eva.evaluation(weightMatrix1,TreeForEachInst1 ,num_vertex2, startVertex2);
			pop.get(i).setFactorialCost(temp);
		}

	}

	public static void sortByCostIndex(ArrayList<Individual> pop, int pop_size, int index) {
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).cost = pop.get(i).getFactorialCost()[index];
		}
		Collections.sort(pop, ChromosomeCmp.compareByFactorialCost); 

	}

	public static void printArray(double[][] arr) {
		for (int i = 0; i < arr.length; i++) {
			for (int j = 0; j < arr[i].length; j++) {
				System.out.print(arr[i][j] + " ");
			}
			System.out.println();
		}

	}
	
	/**
	 * this method to compare scalar fitness of two individuals 
	 * @param x : individual 1
	 * @param y : individual 2
	 * @return 1,0,-1 corresponding the value of scalar fitness 
	 */
	public static int sortByScalarFitness(Individual x, Individual y){
		
		if(x.getScalarFitness() < y.getScalarFitness()){
			return 1;
			
		}else if( x.getScalarFitness() > y.getSkillFactor()){
	
	       return -1;		
		}else{
			return 0;
		}

	}
	
	public static String getDateFromMillis(long millis) {
        String string = String.format("%02d:%02d:%02d.%03d",
                TimeUnit.MILLISECONDS.toHours(millis), TimeUnit.MILLISECONDS.toMinutes(millis) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(millis)),
                TimeUnit.MILLISECONDS.toSeconds(millis) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis)), millis - TimeUnit.SECONDS.toMillis(TimeUnit.MILLISECONDS.toSeconds(millis)));
        return string;
    }

	
	
	
	/**
	 * initialize population 
	 * @param clusters : the set of vertices 
	 * @param popLength  
	 * @param rnd
	 * @return
	 */

	public ArrayList<Individual>  initiaizePopulation(ArrayList<Cluster> clusters, int popLength, Random rnd){
		   ArrayList<Individual> population = new ArrayList<Individual>();
		for(int m = 0; m < popLength; m++){
		ArrayList<Integer> temp = new ArrayList<Integer>();
		ArrayList<Integer> temp1 = new ArrayList<Integer>();
		ArrayList<Integer> temp3 = new ArrayList<Integer>();
		temp3.removeAll(temp3);

		Individual individual = new Individual();
	 
	
	int numberOfClusters = clusters.size();
	int[] presentVertex = new int[numberOfClusters];
	
	for( int  i = 0; i < numberOfClusters; i ++ ){
	
	temp.add(i,i);
	int numberClusterVertex = clusters.get(i).getCluster().size();
	 presentVertex[i] =  rnd.nextInt(InitializeChromosome.minClusterVertices[i]);
	for(int j = 0; j < numberClusterVertex; j ++){
		temp1.add(j,j);
	}
	Collections.shuffle(temp1, rnd);
	temp1 = deletePruferNumber(temp1, rnd);
	temp3.addAll(temp1);
	temp1.removeAll(temp1);
	
	}
	Collections.shuffle(temp, rnd);
	temp = deletePruferNumber(temp, rnd);
	temp3.addAll(temp);
	temp.removeAll(temp);
	int num  = temp3.size() + numberOfClusters;
	int k = 0;
	for(int i = 0; i < num; i ++){ 
		if(i < temp3.size()){
		individual.getGene1()[i] = temp3.get(i);
		} else{
			individual.getGene1()[i] = presentVertex[k];
			k++;
		}
	}
	population.add(individual);
		}
		return population;
	}
	
	
/**
 * delete randomly  from the number  at first code. delete two last elements
 * @param individual : the number at first 
 * @param rnd : random seed
 * @return the prufer number code after delete
 */
	
	public  ArrayList<Integer>  deletePruferNumber(ArrayList<Integer> individual, Random rnd){
		int geneLength = individual.size();
		individual.remove(geneLength - 1);
		individual.remove(geneLength - 2);
		return individual;
	}
	/**
	 * calculate the length of prufernumber gene;
	 * @param maxClusters
	 * @return
	 */
	public  int calCulateGenLength( ArrayList<Cluster> maxClusters){
	int numberOfCluster  = maxClusters.size();
	int geneLength = 0;
	for( int i = 0; i < numberOfCluster; i++){
		int numberClusterVertices =  maxClusters.get(i).getCluster().size();
		geneLength += (numberClusterVertices - 2);
	}
	
	geneLength += (2*numberOfCluster - 2);
	
	return geneLength;
	}
	
	public  double calMuRate( ArrayList<Cluster> maxClusters){
		
		int numberOfCluster  = maxClusters.size();
		int numberVertices = 0;
		for( int i = 0; i < numberOfCluster; i++){
			int numberClusterVertices =  maxClusters.get(i).getCluster().size();
			numberVertices += (numberClusterVertices);
		}
		return 1.0/numberVertices;
		
		
	}
	
	
}
