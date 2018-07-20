/* *********************************************************************
 * ********* Author¡¯s name(s): Yifan Chen
 * Course Title: Artificial Intelligence
 * Semester: Fall 2017
 * Assignment Number 4
 * Submission Date: 11/26/2017
 * Purpose: This program implement the IBeam problem using Genetic Algorithm with importing 
 * 			jenetics-3.9.0.jar API. The 2D diagram is implemented by PlotPackage.jar API.
 * Input: java GA
 * Output: The value of x1, x2, x3 and x4 which belongs to the best fitness in every generation.
 * 		   Four 2D diagrams: 1) Diagram of the best f1 in each generation
 * 							 2) Diagram of the average f1 in the entire generation
 * 							 3) Diagram of the best f2 in each generation
 * 							 4) Diagram of the average f2 in the entire generation
 * Help: I worked alone.
 * ************************************************************************
 * ****** */

import static org.jenetics.engine.EvolutionResult.toBestPhenotype;

import java.awt.Color;
import java.util.Iterator;
import javax.swing.WindowConstants;

import org.jenetics.BitChromosome;
import org.jenetics.BitGene;
import org.jenetics.Genotype;
import org.jenetics.Mutator;
import org.jenetics.Optimize;
import org.jenetics.Phenotype;
import org.jenetics.SinglePointCrossover;
import org.jenetics.engine.Engine;
import org.jenetics.engine.EvolutionResult;

import jahuwaldt.plot.LinearAxisScale;
import jahuwaldt.plot.Plot2D;
import jahuwaldt.plot.PlotAxis;
import jahuwaldt.plot.PlotPanel;
import jahuwaldt.plot.PlotWindow;
import jahuwaldt.plot.SimplePlotXY;

public class GA {
	// the boundary of each gene, x1, x2, x3, x4
    private static final double X1_MIN = 10;
    private static final double X1_MAX = 80;
    private static final double X2_MIN = 10;
    private static final double X2_MAX = 50;
    private static final double X3_MIN = 0.9;
    private static final double X3_MAX = 5.0;
    private static final double X4_MIN = 0.9;
    private static final double X4_MAX = 5.0;
    
    // constants: population size, generation size
    private static final int POP_SIZE = 50;
    private static final int GENERATIONS = 100;
    
    // crossover probability, mutation probability
    private static final double crossoverProb = 0.75;
    private static final double mutationProb = 0.001;
    
    // parameters: alpha, beta
    private static final double evalA = 0;
    private static final double evalB = 1;
    
    // xArr: diagram x axis£¬ yArr: diagram y axis
    private static double[] xArr = new double[GENERATIONS];
	private static double[] yArr = new double[GENERATIONS];
	
	// arrays of average f1, average f2, best f1, best f2 
    private static double avgF1Generation[];
    private static double avgF2Generation[];
    private static double bestF1Generation[];
    private static double bestF2Generation[];
    
    // length of each gene
    private static int x1DigitNum,x2DigitNum,x3DigitNum,x4DigitNum;
    
    // initializes variables
    private static void init() {
        // calculate the length of each gene
    	x1DigitNum = geneDigitNum(X1_MIN, X1_MAX);
    	x2DigitNum = geneDigitNum(X2_MIN, X2_MAX);
    	x3DigitNum = geneDigitNum(X3_MIN, X3_MAX);
    	x4DigitNum = geneDigitNum(X4_MIN, X4_MAX);
    	
    	avgF1Generation = new double[GENERATIONS];
    	avgF2Generation = new double[GENERATIONS];
    	bestF1Generation = new double[GENERATIONS];
    	bestF2Generation = new double[GENERATIONS];
    }
    
    //****************************************************** 
  	//*** Purpose: calculate digit number of each gene
  	//*** Input:  the minimum and maximum of gene
  	//*** Output: the number of bits of each gene
  	//******************************************************
    public static int geneDigitNum(double min, double max){
    	double d=0.0;
    	int m=0;
    	d=(max-min)*10000;
    	m=(int) Math.ceil(Math.log(d)/Math.log(2));
    	return m;
    }
    
 	//****************************************************** 
  	//*** Purpose: change the binary array into number based on ten
  	//*** Input:  the binary array
  	//*** Output: the number based on ten
  	//******************************************************
    private static double binaryArrayToNum(byte[] binaryArray) {  
    	int result = 0;
        for (int i = binaryArray.length-1, k=0; i >=0 ; i--, k++)
            if (binaryArray[i] == 1)
                result += Math.pow(2, k);
        return result;  
    }
    
    // ********************************************************************
    // *** Purpose: decode the BitChromosome, change the BitChromosome into double value
    // *** Input: A BitChromosome, the lower and upper bounds of the gene, and its length
    // *** Output: The corresponding double value
    // ********************************************************************
    private static double decode(BitChromosome bitChromosome, double lowerBound, double upperBound, int bits) {
        if (bitChromosome == null) throw new IllegalArgumentException("bitChromosome is null!");
        if (upperBound < lowerBound) throw new IllegalArgumentException("invalid boundary arguments!");
        if (bits <= 0) throw new IllegalArgumentException("illegal bits!");

        byte[] geneArray=new byte[bits];
        int index = 0;
        double result=0.0;
        Iterator<BitGene> it = bitChromosome.iterator();

        while (it.hasNext()) {
        	BitGene bitGene = it.next();
        	if (bitGene.booleanValue()==true)
        		geneArray[index++]=1;
        	else
        		geneArray[index++]=0;
        }
        
        result=binaryArrayToNum(geneArray);
        result = lowerBound + result * (Math.abs(upperBound - lowerBound) / (Math.pow(2, bits) - 1));
        return result;
    }
    
    // ********************************************************************
    // *** Purpose: according to parameters x1, x2, x3 and x4, calculate the f1 value
    // *** Input: the value of x1, x2, x3 and x4
    // *** Output: the calculated valued of f1
    // ********************************************************************
    public static double f1(double x1,double x2,double x3,double x4){
    	double value=2*x2*x4+x3*(x1-2*x4);
    	return value;
    }
    
    // ********************************************************************
    // *** Purpose: according to parameters x1, x2, x3 and x4, calculate the f2 value
    // *** Input: the value of x1, x2, x3 and x4
    // *** Output: the calculated valued of f2
    // ********************************************************************
    public static double f2(double x1,double x2,double x3,double x4){
    	double value=60000/(x3*(x1-2*x4)*(x1-2*x4)*(x1-2*x4)+2*x2*x4*(4*x4*x4+3*x1*(x1-2*x4)));
    	return value;
    }
    
    // ********************************************************************
    // *** Purpose: according to parameters x1, x2, x3 and x4, calculate the fitness value
    // *** Input: the value of x1, x2, x3 and x4
    // *** Output: the calculated valued of fitness
    // ********************************************************************
    public static double calculateFitness(double x1,double x2,double x3,double x4){
    	return evalA*f1(x1,x2,x3,x4)+evalB*f2(x1,x2,x3,x4);
    }
    
    // calculate the fitness using assigned weights
    private static double fitness(Genotype<BitGene> genotype) {
        double x1, x2, x3, x4, fitness;

        if (genotype == null) throw new IllegalArgumentException("genotype is null!");

        // decode the chromosomes
        x1 = decode(genotype.getChromosome(0).as(BitChromosome.class), X1_MIN, X1_MAX, x1DigitNum);
        x2 = decode(genotype.getChromosome(1).as(BitChromosome.class), X2_MIN, X2_MAX, x2DigitNum);
        x3 = decode(genotype.getChromosome(2).as(BitChromosome.class), X3_MIN, X3_MAX, x3DigitNum);
        x4 = decode(genotype.getChromosome(3).as(BitChromosome.class), X4_MIN, X4_MAX, x4DigitNum);

        // calculate the fitness using assigned weights
        fitness = calculateFitness(x1, x2, x3, x4);
        
        return fitness;
    }
    
    // ********************************************************************
    // *** Purpose: calculate the average f1, f2 value and best f1, f2 in each generation by using plotting API
    // *** Input: the EvolutionResult of type <BitGene, Double>
    // *** Output: Each geometric side constraints and the overall fitness value
    // ********************************************************************
    private static void update(final EvolutionResult<BitGene, Double> result) {
        if (result == null) throw new IllegalArgumentException("result is null!");

        int generation = (int) result.getGeneration();
        int popSize = result.getPopulation().size();
        double x1, x2, x3, x4;
        double bestF1, bestF2, sumF1 = 0.0, sumF2 = 0.0;
        double avg;
        Genotype<BitGene> bestGeno = result.getBestPhenotype().getGenotype();

        // calculate the best geno's f1 and f2
        x1 = decode(bestGeno.getChromosome(0).as(BitChromosome.class), X1_MIN, X1_MAX, x1DigitNum);
        x2 = decode(bestGeno.getChromosome(1).as(BitChromosome.class), X2_MIN, X2_MAX, x2DigitNum);
        x3 = decode(bestGeno.getChromosome(2).as(BitChromosome.class), X3_MIN, X3_MAX, x3DigitNum);
        x4 = decode(bestGeno.getChromosome(3).as(BitChromosome.class), X4_MIN, X4_MAX, x4DigitNum);

        bestF1 = 2*x2*x4+x3*(x1-2*x4);
        bestF2 = 60000/(x3*(x1-2*x4)*(x1-2*x4)*(x1-2*x4)+2*x2*x4*(4*x4*x4+3*x1*(x1-2*x4)));
        
        // test
//        System.out.printf("f1 = %.6f f2=%.6f", bestF1,bestF2);
//        System.out.println();

        bestF1Generation[generation - 1] = bestF1;
        bestF2Generation[generation - 1] = bestF2;

        // print the best fitness info to the console
        System.out.println("Generation " + generation + "'s Best:");
        System.out.println("-------------------------------------------------------------------------------");
        System.out.printf("x1 = %.6f\tx2 = %.6f\tx3 = %.6f\tx4 = %.6f\n", x1, x2, x3, x4);
        System.out.printf("Fitness: %.6f\n", result.getBestFitness());
        System.out.println();

        Iterator<Phenotype<BitGene, Double>> it = result.getPopulation().iterator();
        while (it.hasNext()) {
            Genotype<BitGene> current = it.next().getGenotype();
            // calculate the average f1 and f2
            x1 = decode(current.getChromosome(0).as(BitChromosome.class), X1_MIN, X1_MAX, x1DigitNum);
            x2 = decode(current.getChromosome(1).as(BitChromosome.class), X2_MIN, X2_MAX, x2DigitNum);
            x3 = decode(current.getChromosome(2).as(BitChromosome.class), X3_MIN, X3_MAX, x3DigitNum);
            x4 = decode(current.getChromosome(3).as(BitChromosome.class), X4_MIN, X4_MAX, x4DigitNum);

            sumF1 += f1(x1, x2, x3, x4);
            sumF2 += f2(x1, x2, x3, x4);
        }

        avg = sumF1 / popSize;
        avgF1Generation[generation - 1] = avg;

        avg = sumF2 / popSize;
        avgF2Generation[generation - 1] = avg;
    }
    
    // ********************************************************************
    // *** Purpose: generate the plotting diagram
    // *** Input: double[] yArr, String title, String xAxisTitle, String yAxisTitle, int locX, int locY
    // *** Output: none
    // ********************************************************************
    public static void createPlot(double[] yArr, String title, String xAxisTitle, String yAxisTitle, int locX, int locY){
		for (int i=0; i<GENERATIONS; i++){
			xArr[i] = i;
			yArr[i] = yArr[i]; // your GASolver should fill up this array
		}
		Plot2D aPlot = new SimplePlotXY(xArr, yArr, title, xAxisTitle, yAxisTitle, null, null, null);
		// Make the horizontal axis a log axis.
		PlotAxis xAxis = aPlot.getHorizontalAxis();
		xAxis.setScale(new LinearAxisScale());
		PlotPanel panel = new PlotPanel(aPlot);
		panel.setBackground( Color.white );
		PlotWindow window = new PlotWindow("SimplePlotXY Plot Window", panel);
		window.setSize(500, 300);
		window.setLocation(locX,locY); // location on screen
		window.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		window.show();
	}
    
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		// initial the variables
		init();

        // define the genotype of four BitChromosomes
        // each chromosome has its calculated length
        final Genotype<BitGene> genotype = Genotype.of(
                BitChromosome.of(x1DigitNum),
                BitChromosome.of(x2DigitNum),
                BitChromosome.of(x3DigitNum),
                BitChromosome.of(x4DigitNum)
        );
        
        // define and configure the engine to fit the requirements
        // of this particular question
        final Engine<BitGene, Double> engine = Engine.builder(
                GA::fitness, genotype)
                .populationSize(POP_SIZE)
                .optimize(Optimize.MINIMUM)
                .alterers(
                        new SinglePointCrossover<>(crossoverProb),
                        new Mutator<>(mutationProb)
                ).build();
        
        // start the engine and run the Genetic algorithm
        engine.stream()
                .peek(GA::update)
                .limit(GENERATIONS)
                .collect(toBestPhenotype());
        
        // show the plotting diagrams
        createPlot(bestF1Generation,"Cross Section Area","Generation No.", "Best Particle Fitness",150,150);
		createPlot(avgF1Generation,"Cross Section Area", "Generation No.", "Avg. Population Fitness",680,150);
		createPlot(bestF2Generation,"Static Deflection", "Generation No.", "Best Particle Fitness",150,480);
		createPlot(avgF2Generation,"Cross Section Area", "Generation No.", "Avg. Population Fitness",680,480);
	}
}
