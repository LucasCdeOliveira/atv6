package problems.qbf.solvers;

import metaheuristics.ga.AbstractGA;
import problems.qbf.KQBF;
import problems.qbf.QBF;
import solutions.Solution;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Metaheuristic GA (Genetic Algorithm) for
 * obtaining an optimal solution to a QBF (Quadractive Binary Function --
 * {@link #QuadracticBinaryFunction}).
 *
 * @author ccavellucci, fusberti
 */
public class GA_KQBF extends AbstractGA<Integer, Integer> {

    /**
     * Constructor for the GA_QBF class. The QBF objective function is passed as
     * argument for the superclass constructor.
     *
     * @param generations  Maximum number of generations.
     * @param popSize      Size of the population.
     * @param mutationRate The mutation rate.
     * @param filename     Name of the file for which the objective function parameters
     *                     should be read.
     * @throws IOException Necessary for I/O operations.
     */
    public GA_KQBF(Integer generations, Integer popSize, Double mutationRate, String filename) throws IOException {
        super(new KQBF(filename), generations, popSize, mutationRate);
    }

    /**
     * {@inheritDoc}
     * <p>
     * This createEmptySol instantiates an empty solution and it attributes a
     * zero cost, since it is known that a QBF solution with all variables set
     * to zero has also zero cost.
     */
    @Override
    public Solution<Integer> createEmptySol() {
        Solution<Integer> sol = new Solution<Integer>();
        sol.cost = 0.0;
        return sol;
    }

    /*
     * (non-Javadoc)
     *
     * @see metaheuristics.ga.AbstractGA#decode(metaheuristics.ga.AbstractGA.
     * Chromosome)
     */
    @Override
    protected Solution<Integer> decode(Chromosome chromosome) {

        Solution<Integer> solution = createEmptySol();
        for (int locus = 0; locus < chromosome.size(); locus++) {
            if (chromosome.get(locus) == 1) {
                solution.add(new Integer(locus));
            }
        }

        ObjFunction.evaluate(solution);
        return solution;
    }

    /*
     * (non-Javadoc)
     *
     * @see metaheuristics.ga.AbstractGA#generateRandomChromosome()
     */
    @Override
    protected Chromosome generateRandomChromosome() {

        Chromosome chromosome = new Chromosome();
        ArrayList<Integer> identityList = new ArrayList<>();

        //Setting identityList to identify the genes from 0 to 99
        for (int i = 0; i < chromosomeSize; i++) {
            identityList.add(i);
        }

        //Shuffling to randomize the order that the genes will be selected
        Collections.shuffle(identityList);

        //Setting all genes to 0
        for (int i = 0; i < chromosomeSize; i++) {
            chromosome.add(0);
        }

        for (int i = 0; i < chromosomeSize; i++) {
            //I only randomize the possible insertion if the gene fits in the chromosome
            if (ObjFunction.shouldInsert(i, decode(chromosome))) {
                chromosome.set(i, rng.nextInt(2));
            } else {
                //otherwise I can't add the gene to the chromossome
                chromosome.set(i, 0);
            }
        }
        return chromosome;

    }

    /*
     * (non-Javadoc)
     *
     * @see metaheuristics.ga.AbstractGA#fitness(metaheuristics.ga.AbstractGA.
     * Chromosome)
     */
    @Override
    protected Double fitness(Chromosome chromosome) {

        return decode(chromosome).cost;

    }

    /*
     * (non-Javadoc)
     *
     * @see
     * metaheuristics.ga.AbstractGA#mutateGene(metaheuristics.ga.AbstractGA.
     * Chromosome, java.lang.Integer)
     */
    @Override
    protected void mutateGene(Chromosome chromosome, Integer locus) {
        Integer oldGene = chromosome.get(locus);
        chromosome.set(locus, 1 - chromosome.get(locus));
        if (ObjFunction.evaluateWeight(decode(chromosome)) <= ((KQBF) ObjFunction).getCapacity())
            //This is a valid mutation because it doesn't exceed the weight
            //So I can change in the original chromosome
            chromosome.set(locus, 1 - chromosome.get(locus));
        else {
            //Otherwise I revert the mutation do nothing because it is not a valid mutation
            chromosome.set(locus, oldGene);
        }
    }

    @Override
    protected Population crossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome parent1 = parents.get(i);
            Chromosome parent2 = parents.get(i + 1);

            int crosspoint1 = rng.nextInt(chromosomeSize + 1);
            int crosspoint2 = crosspoint1 + rng.nextInt((chromosomeSize + 1) - crosspoint1);

            Chromosome offspring1 = new Chromosome();
            Chromosome offspring2 = new Chromosome();

            for (int j = 0; j < chromosomeSize; j++) {
                if (j >= crosspoint1 && j < crosspoint2) {
                    offspring1.add(parent2.get(j));
                    offspring2.add(parent1.get(j));
                } else {
                    offspring1.add(parent1.get(j));
                    offspring2.add(parent2.get(j));
                }
            }

            //If the parents can't generate childs with valid weight I just copy the parents to the next population
            //Since I know the parents have valid weights

            if (ObjFunction.evaluateWeight(decode(offspring1)) <= ((KQBF) ObjFunction).getCapacity())
                //Child 1 have a valid weight
                offsprings.add(offspring1);
            else
                //Child 1 doesn't have a valid weight
                offsprings.add(parent1);

            if (ObjFunction.evaluateWeight(decode(offspring2)) <= ((KQBF) ObjFunction).getCapacity())
                //Child 2 have a valid weight
                offsprings.add(offspring2);
            else
                //Child 2 doesn't have a valid weight
                offsprings.add(parent2);
        }

        return offsprings;
    }

    //remove elements until the chromosome be valid
    private Chromosome validateChromosome(Chromosome chromosome){
        while(ObjFunction.evaluateWeight(decode(chromosome)) > ObjFunction.getCapacity()){
            chromosome.set(getWorstGene(chromosome), 0);
        }
        return chromosome;
    }

    //get the worst gene from chromosome
    private Integer getWorstGene(Chromosome chromosome){
        Solution sol = decode(chromosome);
        Double initialWeight = ObjFunction.evaluateWeight(sol);
        Integer worstGene = 0;
        Double worstRate = Double.MAX_VALUE;
        for (int i = 0; i < chromosomeSize; i++) {
            Chromosome tempChromosome = new Chromosome();
            tempChromosome.addAll(chromosome);
            //removing element i
            tempChromosome.set(i, 0);
            Solution newSol = decode(tempChromosome);
            Double newWeight = ObjFunction.evaluateWeight(newSol);
            Double newRate = newSol.cost/newWeight;
            if(newRate < worstRate){
                worstRate = newRate;
                worstGene = i;
            }
        }
        return worstGene;
    }

    //Here instead of useing 2 crosspoints we will randomly decide wich genes will come from parent1 and parent2
    protected Population universalCrossover(Population parents) {

        Population offsprings = new Population();

        //if 0 we will get from parent1, if 2 we will get from parent2 for offspring1
        //if 0 we will get from parent2, if 2 we will get from parent1 for offsprint2
        ArrayList<Integer> pattern = new ArrayList<>();
        for (int i = 0; i < chromosomeSize; i++){
            pattern.add(rng.nextInt(2));
        }

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome parent1 = parents.get(i);
            Chromosome parent2 = parents.get(i + 1);

            Chromosome offspring1 = new Chromosome();
            Chromosome offspring2 = new Chromosome();

            for (int j = 0; j < chromosomeSize; j++) {
                if (pattern.get(j) == 0) {
                    offspring1.add(parent1.get(j));
                    offspring2.add(parent2.get(j));
                } else {
                    offspring1.add(parent2.get(j));
                    offspring2.add(parent1.get(j));
                }
            }

            offsprings.add(validateChromosome(offspring1));
            offsprings.add(validateChromosome(offspring2));

        }

        return offsprings;

    }

    //Because of I typed the chromossome to Integer to validate it
    @Override
    public Solution<Integer> solve() {

        /* starts the initial population */
        Population population = initializePopulation();

        bestChromosome = getBestChromosome(population);
        bestSol = decode(bestChromosome);
        System.out.println("(Gen. " + 0 + ") BestSol = " + bestSol);

        /*
         * enters the main loop and repeats until a given number of generations
         */
        for (int g = 1; g <= generations; g++) {

            Population parents = selectParents(population);

            Population offsprings = universalCrossover(parents);

            Population mutants = mutate(offsprings);

            Population newpopulation = selectPopulation(mutants);

            population = newpopulation;

            bestChromosome = getBestChromosome(population);

            if (fitness(bestChromosome) > bestSol.cost) {
                bestSol = decode(bestChromosome);
                if (verbose)
                    System.out.println("(Gen. " + g + ") BestSol = " + bestSol);
            }

        }

        return bestSol;
    }

    /**
     * A main method used for testing the GA metaheuristic.
     */
    public static void main(String[] args) throws IOException {

        long startTime = System.currentTimeMillis();
        GA_KQBF ga = new GA_KQBF(1000, 100, 1.0 / 100.0, "instances/kqbf/kqbf200");
        Solution<Integer> bestSol = ga.solve();
        System.out.println("maxVal = " + bestSol);
        KQBF evaluateCost = new KQBF("instances/kqbf/kqbf200");
        System.out.println("weight of solution = " + evaluateCost.evaluateWeight(bestSol));
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

    }

}
