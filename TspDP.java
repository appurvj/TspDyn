import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.math.BigDecimal;
public class TspDP{
	public static final int DOUBLE_SCALE = 7;
	public static final int EARTH_RADIUS = 6371; // source Wikipedia
	private double [][] dist;
	private int N;
	private double C[][];
	private int[] optPath;


	public static double[][] readFile(String path){
    ArrayList<double[]> coordArrList;
		try{
			Scanner scan  = new Scanner(new File(path));
			coordArrList = new ArrayList<double[]>();
			Pattern xCoord = Pattern.compile("\\(\\s*(-*\\d*\\.\\d*)\\s*,\\s*");
			Pattern yCoord = Pattern.compile(",\\s*(-*\\d*\\.\\d*)\\s*\\)");

			while(scan.hasNextLine()){
				String line;
				if((line = scan.nextLine() )== "")
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
//System.out.println(coord[0] + ",  " + coord[1]);
			}

		}catch (FileNotFoundException e){
			System.out.println("File Not Found. Kindly check file path provided");
			return null;
		}
		int nVert = coordArrList.size();
		double[][] distMat = new double[nVert+1][nVert+1];
		for(int i = 1; i < nVert; i++){
			for(int j = i+1; j <= nVert; j++){
				double[] c1 = coordArrList.get(i-1);
				double[] c2 = coordArrList.get(j-1);
//				double[] diff = {c2[0]-c1[0], c2[1] - c1[1]};
				distMat[i][j] = distMat[j][i] = haversineDist(c1,c2);
					//Math.sqrt(diff[0]*diff[0] + diff[1]*diff[1]); //(new BigDecimal(Math.sqrt(diff[0]*diff[0] + diff[1]*diff[1]))).setScale(DOUBLE_SCALE, BigDecimal.ROUND_HALF_UP).doubleValue();
			}
		}
/*
for(int i = 0; i <= nVert; i++){
	for(int j = 0; j <=nVert; j++)
		System.out.print(distMat[i][j] + "    ");
	System.out.println();
}
*/
		return distMat;

	}

	public static double haversineDist(double[] c1, double[] c2 ){
		double latDist = (c2[0] - c1[0]) * Math.PI/180;
		double longDist = (c2[1] - c1[1]) * Math.PI/180;
		double a = Math.sin(latDist/2.0) * Math.sin(latDist/2.0) + 
							 Math.cos(c1[0] * Math.PI/180)*Math.cos(c2[0]*Math.PI/180)*
							 Math.sin(longDist/2.0) * Math.sin(longDist/2.0);
		double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
		 return c * EARTH_RADIUS;


	}

	public TspDP(double[][] distMat){ //size of distMat is (N+1, N+1)
		this.dist = distMat;
		this.N = distMat.length - 1;
		this.C = new double[N+1][1<<(N-1)];
		this.populateCostMatrix();
		this.calcOptPath();
	}

	private void populateCostMatrix(){
		for(int i  = 2; i <= N; i++){
			C[i][(1<<(i-2))] = dist[1][i];
		}
		
		for(int nV = 2; nV < N ; nV++){ //number of vertices in perm (apart from 1)
			for(int s = (1<<(nV))-1; s < 1<<(N-1); s = nextNumSameBits(s)){ //perm for vertices selected
				for(int k = 2; k <= N; k++){ //ending vertex
					if((s&(1<<(k-2))) != 0){ // proceed with calculation only if the required vertex is part of the perm
					  C[k][s] = Double.MAX_VALUE;
					  for(int l = 2; l <= N; l++){ //go through various routing possibilities ending in vertex k for perm s
					  	if(k!=l && (s&(1<<(l-2)))!= 0){
			 		  		C[k][s] = Math.min(C[l][s^(1<<(k-2))]+dist[l][k], C[k][s]);	
							}
						}
					}
				}
			}
		}

/*
System.out.println("Printing C matrix");
for(int k = 2; k <=N; k++){
	System.out.println("k = " + k);
	for(int j = 0; j < 1<<(N-1); j++)
		System.out.print(C[k][j] + "     ");
	System.out.println("\n\n");
}
*/


	}

	public void calcOptPath(){
		int vIdx = N-1;//we will calculate the path from the last vertex
		this.optPath = new int[N]; 
		optPath[0] = 1;//since the path will always start from the firstVertex;
		int currS = (1<<(N-1))-1, prevV = 0;//Since dist[*][0] = 0
		while(vIdx > 0){
			double currMin = Double.MAX_VALUE;
		  for(int k = 2; k <=N; k++){
				if(((currS & (1<<(k-2))) !=0) && currMin > C[k][currS]+dist[k][prevV]){
				  currMin = C[k][currS] + dist[k][prevV];
					optPath[vIdx] = k;	
				}	
		  }
			prevV = optPath[vIdx];
			currS^=1<<(prevV-2);
			vIdx--;
		}

	}

	public int[] getOptPath(){
		return optPath;
	}

  public void printOptPath(){
    for(int i = 0; i < N; i++)
	    System.out.println(optPath[i]);
	}
	
	public static int nextNumSameBits(int x){
		return (((x^(x + (x&(-x))))/(x&(-x)))>>2)|(x + (x&(-x)));
	}

	public static void Usage(){
		System.out.println("Usage: java TspDP <path to file>");
	}

	public static void main(String[] args){
		/*
		double[][] d = {{0.0,0.0,0.0,0.0,0.0},{0.0, 0.0, 20.0, 42.0, 35.0},
		 {0.0, 20.0, 0.0, 30.0, 34.0},
		 {0.0, 42.0, 30.0, 0.0, 12.0},
		 {0.0, 35.0, 34.0, 12.0, 0.0}};

		TspDP tsp = new TspDP(d);
		tsp.calcOptPath();
    */
		if(args.length == 1){
    	double[][] d = readFile(args[0]);
			if(d != null){
				TspDP tsp = new TspDP(d);
				tsp.printOptPath();
			}
		}else
			TspDP.Usage();
	}
}
