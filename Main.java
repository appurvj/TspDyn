import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.math.BigDecimal;

/**
 * 							Main.java
 * This class provides the functionality of, given vertex and distances between 
 * them, computing the optimal tour starting from the first vertex and ending at
 * some other vertex. One must note that this is a slightly modified version of 
 * the original Traveling Salesman problem because the tour does not go back to 
 * the starting vertex. This may lead to different paths than traditional TSP
 * 
 * The approach used to solve the problem involves using Dynamic Programming 
 * and populating the cost matrix and then back tracing the optimal path using
 * the populated cost matrix
 * 
 * Note: Since integers are used asbitsets have been used to represent path 
 * permutations, this program can only work for problems upto 30 - 31 vertices
 * 
 * @author Appurv Jain
 * @sources http://www.lsi.upc.edu/~mjserna/docencia/algofib/P07/dynprog.pdf
 */

public class Main{
  public static final int EARTH_RADIUS = 6371; // source https://www.princeton.edu/~achaney/tmve/wiki100k/docs/Earth_radius.html
  public static final int MAX_VERTICES = 30;// keeping it to 30 to be safe since using integer bitsets
  private final double [][] dist; //Since vertices are GPS coordinates, haversine distance is used
  private final int N;
  private double C[][]; //cost matrix
  private final int[] optPath; //Optimal minimum cost path generated will be stored here


 
  public int[] getOptPath(){
    return optPath;
  }

  
  public void printOptPath(){
    for(int i = 0; i < N; i++)
      System.out.println(optPath[i]);
  }
  
  /**
   * Constructor accept the the distance matrix and also populates the cost
   * matrix and computes the optimal path. The calculation is done in the
   * constructor itself because this object is useful only after that calculation
   * @param distMat : It is imp to note that for the sake of this problem,
   * the distMat is assumed to have and extra row and column for ease of dereferencing
   */
  public Main(double[][] distMat){ //size of distMat is (N+1, N+1)
    this.dist = distMat;
    this.N = distMat.length - 1;
    this.C = new double[N][1<<(N-1)];
    this.populateCostMatrix();
    this.optPath = this.calcOptPath();
  }
  
  
  /**
   * This is one of the main methods where most of the computation happens.
   * 
   * Every path begins from 1
   * C[k][s] represents and path from 1 ending in vertex k, using the vertices
   * represent by the bitset s
   * 
   * One of the optimizations in this method is the way we progress through the
   * various permutations.
   * 
   * We progress through permutations first involving 2 vertices then 3 and so
   * on. The nextumSameBits method helps us progress to the next sequence 
   * involving nV vertices in constant time using bit manipulation. 
   * This helps us save unnecessary wasteful looping
   * 
   * To populate a cell in the matrix, one needs all the cells corresponding to
   * a tour with lesser number of vertices populated
   */
  private void populateCostMatrix(){
    for(int i  = 2; i <= N; i++)
      C[i-1][(1<<(i-2))] = dist[1][i];
    
    int sNew;
    double tempDist;
    for(int nV = 2; nV < N ; nV++){ //number of vertices in perm (apart from 1)
      for(int s = (1<<(nV))-1; s < 1<<(N-1); s = nextNumSameBits(s)){ //permutation for vertices selected
        for(int k = 2; k <= N; k++){ //ending vertex
          if((s&(1<<(k-2))) != 0){ // proceed with calculation only if the required vertex is part of the perm
            C[k-1][s] = Double.MAX_VALUE;
            //we need to remove k from the perm and search previously calculated paths
            sNew = s^(1<<(k-2)); 
            for(int l = 2; l <= N; l++) //go through various routing possibilities ending in vertex k for perm s
              if((sNew&(1<<(l-2)))!= 0 && (tempDist = C[l-1][sNew]+dist[l][k]) < C[k-1][s])
                 C[k-1][s] = tempDist;  
          }
        }
      }
    }
  }
  
  
  /**
   * This method uses the populates cost matrix to calculate the optimal path
   * @return : array containing the order in which vertices should be visited
   */
  private int[] calcOptPath(){ 
    int[] localOptPath = new int[N];
    localOptPath[0] = 1;//since the path will always start from the firstVertex;
    int vIdx = N-1;//we will calculate the path backwards
    int currS = (1<<(N-1))-1, prevV = 0;//Starting point since, since dist[*][0] = 0
    double currMin, temp;
    while(vIdx > 0){
      currMin = Double.MAX_VALUE;
      for(int k = 2; k <=N; k++)
        if(((currS & (1<<(k-2))) !=0) && (temp = C[k-1][currS]+dist[k][prevV]) < currMin){
          currMin = temp;
          localOptPath[vIdx] = k; //curr optimal pah ending in k 
        }  
      prevV = localOptPath[vIdx];
      currS^=1<<(prevV-2);
      vIdx--;
    }
    return localOptPath;
  } 
 
  /**
   * This method is used to parse the file provided and generate a distance
   * matrix. Method is static since independent of class variables
   * @param path: String providing path to the file
   * @return : Distance matrix with an extra row and column
   */
  public static double[][] readFile(String path){
    Scanner scan;
		try{
      scan  = new Scanner(new File(path));
    }catch (FileNotFoundException e){
      System.out.println("File Not Found. Kindly check file path provided");
      return null;
    }
    ArrayList<double[]> coordArrList = new ArrayList<double[]>();
    Pattern xCoord = Pattern.compile("\\(\\s*(-*\\d*\\.\\d*)\\s*,\\s*");
    Pattern yCoord = Pattern.compile(",\\s*(-*\\d*\\.\\d*)\\s*\\)");

    String line;
    while(scan.hasNextLine()){
      line = scan.nextLine();
      if(line.trim().equals(""))
        continue;
      Matcher mX = xCoord.matcher(line);
      Matcher mY = yCoord.matcher(line);
      if(!mX.find() || !mY.find()){
        System.out.println("Kindly ensure the text in the provided file is in the correct format");
        return null;
      }
      double[] coord = new double[2];
      coord[0] = (new BigDecimal(mX.group(1))).doubleValue();
        coord[1] = (new BigDecimal(mY.group(1))).doubleValue();
      coordArrList.add(coord);
    }
    
    int nVert = coordArrList.size();
    if(nVert > MAX_VERTICES){
      System.out.println("Program cannot accept probelms involving more than 30 vertices");
      return null;
    }
    double[][] distMat = new double[nVert+1][nVert+1];
    for(int i = 1; i < nVert; i++){
      for(int j = i+1; j <= nVert; j++){
        double[] c1 = coordArrList.get(i-1);
        double[] c2 = coordArrList.get(j-1);
        distMat[i][j] = distMat[j][i] = haversineDist(c1,c2);
      }
    }
    return distMat;
  }

  /**
   * This method computes the harversine distance between two given GPS coordinates.
   * This particular distance is used because the locations of the offices have
   * provided in latitude and longitude. Method is static because it does not 
   * depend on class variables
   * @param c1: coordinate 1
   * @param c2: coordinate 2
   * @return : Distance in km
   */
  public static double haversineDist(double[] c1, double[] c2 ){
    double latDist = (c2[0] - c1[0]) * Math.PI/180;
    double longDist = (c2[1] - c1[1]) * Math.PI/180;
    double a = Math.sin(latDist/2.0) * Math.sin(latDist/2.0) + 
               Math.cos(c1[0] * Math.PI/180)*Math.cos(c2[0]*Math.PI/180)*
               Math.sin(longDist/2.0) * Math.sin(longDist/2.0);
    double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    return c * EARTH_RADIUS;
  }


  /**
   * This method accepts and integer and returns a higher integer with the same
   * number of bits set.
   * 
   * Therefore for 1 (0001), it will return 2 (0010)
   * 
   * Once again, method is static because it is independent
   * of class variables
   * @param x: Input integer
   */
  private static int nextNumSameBits(int x){
		int rightMostOne = x&(-x); //gives a number with jsut righmost bit set
		int setHigherBit = x + rightMostOne; //set the next higher bit to the rightmost bit
	  int adjRightMostOnesChain = ((x^setHigherBit)/rightMostOne)>>2;
		return setHigherBit | adjRightMostOnesChain;
  }

  private static void Usage(){
    System.out.println("Usage: java Main <path to file>");
  }

  public static void main(String[] args){
    if(args.length == 1){
      double[][] d = readFile(args[0]);
      if(d != null){
      	//long t = System.currentTimeMillis();
        Main tsp = new Main(d);
        tsp.printOptPath();
	      //System.out.println("Time taken: " + (System.currentTimeMillis() - t) + " msec");
      }
    }else
      Main.Usage();
  }
}

