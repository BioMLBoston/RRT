package RRT;

import io.DCDIO;
import io.plainLogger;

import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.logging.Logger;
import pdb.PDBAtom ;
import pdb.PDBMolecule;

import RRT.BBMolecule;
/**
 * @author Dong Luo 2012-07-02
 *
 */


public class RRT {
	private int angleStep=5,  scoreThreshold=600;
        private double closeDistance=0.8;
      	private double acceptConstant = 0.01;
	private RRTNode sNode=null, tNode=null;
	private ArrayList<RRTNode> sTree = new ArrayList<RRTNode>(); private int maxSampleSize = 20;
	private ArrayList<RRTNode> treePool = new ArrayList<RRTNode>();
	private ArrayList<RRTNode> tTree = new ArrayList<RRTNode>();
	private ArrayList<RRTNode> rTree = new ArrayList<RRTNode>();
	ArrayList<Double> minDists = new ArrayList<Double>(), maxDists = new ArrayList<Double>();
	
	private RRTNode sLinkNode = null, tLinkNode = null;
	private String log = "", base = null;
	Logger logger = null;
	DataOutputStream out = null;
	
        
        //Extract the element of BBMolecule from source and target molecules
	public RRT(BBMolecule source, BBMolecule target, String base, int scoreThreshold) {
		this.scoreThreshold = scoreThreshold;
		BBMolecule.setElementsAndFreeAngles(source, target, angleStep);
/*
		String info = "Rigid elments size: "+target.getElementList().size();
		info += "; Free angle number: "+target.getFreeAngleList().size();
		info += "; Free dihedral angle number: "+target.getFreeDihedralAngleList().size();
		System.out.println(info);
		*/
		sNode = new RRTNode(source, source.getRMSD(target));
		tNode = new RRTNode(target, target.getRMSD(target));
		sTree.add(sNode); tTree.add(tNode); treePool.add(sNode);
		base += "/";
		this.base = base;
		System.out.println(String.format("Score threshold: %.2f; Initial score: %.2f\n",
					(float)scoreThreshold, source.getRMSD(target)));
		ArrayList<PDBAtom> sal = source.getPDBAtomList(), tal = target.getPDBAtomList();
		for (int i=0;i<sal.size()-1;i++) {
			double sd = sal.get(i).distance(sal.get(i+1));
			double td = tal.get(i).distance(tal.get(i+1));
			minDists.add(java.lang.Math.min(sd, td));
			maxDists.add(java.lang.Math.max(sd, td));
		}
	}
	
	private RRTNode nearestNode(ArrayList<RRTNode> tree, RRTNode node) {
		double min = 0;
		RRTNode ret = null;
		for (RRTNode n:tree) {
			double d = getSquareDistance(n, node);
			if (d < 0) continue;
			if (ret == null || d < min) {
				min = d;
				ret = n;
			}
		}
		return ret;
	}
	
	private double nearestDistance(ArrayList<RRTNode> tree, double sv) {
		double min = 0;
		RRTNode ret = null;
		for (RRTNode n:tree) {
			double d = Math.abs(n.score- sv);
			if (ret == null || d < min) {
				min = d;
				ret = n;
			}
		}
		return java.lang.Math.sqrt(min);
	}

	private double getDelta(int index, Double dist) {
		double minD = minDists.get(index), maxD = maxDists.get(index);
		if (dist<minD) return minD-dist;
		if (dist>maxD) return maxD-dist;
		return 0;
	}

	private void correctAtomDistance(BBMolecule mol) {
		// correct any distance error caused by numeric error in rotating calculation
		ArrayList<PDBAtom> l = mol.getPDBAtomList();
		for (int j=0;j<l.size()-1;j++) {
			double delta = getDelta(j, l.get(j).distance(l.get(j+1)));
			if (delta<-0.2||delta>0.2) 
                            mol.translate(j, delta);
		}
	}
	
	private RRTNode extendRRT(ArrayList<RRTNode> tree, RRTNode to, boolean bias) {
            
		double steps = 36;
		int maxNotImprovedSteps = 2;
		ArrayList<Double> diff_vector = new ArrayList<Double>();
		ArrayList<Double> step_vector = new ArrayList<Double>();
		RRTNode from = nearestNode(tree, to), min = null;
		BBMolecule m1 = from.getMolecule(), m2 = to.getMolecule(), target = tNode.getMolecule();
                
		double parentScore = from.score, newScore = 0, minScore = -1;
		//ArrayList<Integer>al = m1.getFreeAngleList();
		ArrayList<Integer>dal = m1.getDihedralAngleList();
		//int s1 = al.size();
                int s2 = dal.size();
		double maxDiff = 0;
		/*for (int i:al) {
			double da = (m2.getAngle(i)-m1.getAngle(i))*180/java.lang.Math.PI;
			if (da>180) da -= 360;
			else if (da<-180) da += 360;
			diff_vector.add(da);
			da = java.lang.Math.abs(da);
			if (maxDiff < da) maxDiff = da;
		}*/
		for (int i:dal) {
			double da = (m2.getDihedralAngle(i)-m1.getDihedralAngle(i))*180/java.lang.Math.PI;
			if (da>180) da -= 360;
			else if (da<-180) da += 360;
			diff_vector.add(-da);	// dihedral angle rotate direction is reversed somehow
			da = java.lang.Math.abs(da);
			if (maxDiff < da) maxDiff = da;
		}
		double a = angleStep/maxDiff;
		//System.out.print(String.format("%.2f %.2f, ", a, maxDiff));
		for (double d:diff_vector) step_vector.add(a*d);
		int notImprovedSteps = 0;
		for (int i=0;i<steps;i++) {
			BBMolecule mol = BBMolecule.newBBMolecule(m1);
		
                        
			
			for (int j=0;j<s2;j++) {
				double d = step_vector.get(j)*(i+1);
				if (d>0.2) {
					mol.rotateAroundBond(dal.get(j), dal.get(j)+1, d);
				}
			}
			
			correctAtomDistance(mol);
                        PDBMolecule PDBMol=mol;
			if (mol.isEnergyValid(PDBMol)) {
			
                                newScore = mol.getRMSD(target);
				if (bias && newScore>parentScore) {
					notImprovedSteps++;
					if (notImprovedSteps>maxNotImprovedSteps) break;
				}
				RRTNode n = new RRTNode(from, newScore, step_vector, i+1);
				tree.add(n);
				if (newScore<minScore || minScore<0) {
					minScore = newScore; min = n;
				}
			} else break;
		}
                
		return min;
	}

	private BBMolecule newRandomMolecule(BBMolecule sMol) {
		if (sMol == null) return null;
		BBMolecule tMol = tNode.getMolecule();
                
		ArrayList<Integer>dal = sMol.getDihedralAngleList();
		int s2 = dal.size();
		ArrayList<Double>deltaAngle = new ArrayList<Double>();
		double sum = 0;
	
                
		for (int i:dal) {
			double d = java.lang.Math.abs(tMol.getDihedralAngle(i)-sMol.getDihedralAngle(i));
			if (d>java.lang.Math.PI) d = java.lang.Math.PI*2 - d;
			sum += d;
			deltaAngle.add(sum);
		}
		BBMolecule molecule = null;
		Random g0=new Random(), g1=new Random(), g2=new Random();
		ArrayList<Integer> toBeSampled = new ArrayList<Integer>();
                PDBMolecule PDBMole= molecule;
		do {
			molecule = BBMolecule.newBBMolecule(sMol);
			//	pick angles to be sampled depending on their difference values stored in deltaAngle
			toBeSampled.clear();
			for (int i=0;i<s2;i++) {
				double d = g0.nextDouble()*deltaAngle.get(s2-1);
				int index = 0;
				for (;index<(s2);index++)
					if (d < deltaAngle.get(index)) break;
				if (!toBeSampled.contains(index)) toBeSampled.add(index);
			}
			for (int i:toBeSampled) {
				
					molecule.rotateAroundBond(dal.get(i), dal.get(i)+1, g2.nextDouble()*angleStep*2-angleStep);
			}
			correctAtomDistance(molecule);
                        PDBMole= molecule;
		} while (!(molecule.isEnergyValid(PDBMole)));
		  
                return molecule;
              
	}

	public ArrayList<BBMolecule> searchPathFromSource() {
            
            // write log file
		logger = new plainLogger().getLogger("cg.rrt", base+angleStep+"."+scoreThreshold+".log");
		String dcdFile = angleStep+"."+scoreThreshold+".dcd";
		out = DCDIO.getOutHandle(base+dcdFile);
		BBMolecule source = sNode.getMolecule(), target = tNode.getMolecule();
		DCDIO.writeHeader(out, source, dcdFile);
		DCDIO.writeData(out, source);
		//double minScore = source.getScore(target);
		double minScore = source.getRMSD(target);
                
                logger.info(String.format("%.2f", minScore ));
                
		int maxNodeNum = 1000, maxFailIteration = 100;
		int failNum = 0;
		boolean bias = false;
		Random g1 = new Random(), g2 = new Random(), g3 = new Random();
		int iter = 0;
		while (sTree.size()<maxNodeNum && failNum<maxFailIteration) {
			bias = (g3.nextInt(3) >0);
			
			// create a new random molecule
			RRTNode pNode = null;
			BBMolecule mol = null;
			double sv=1000 ;
			int rIter = 0, sIter = 0;
			while (rIter < maxFailIteration && sIter < maxFailIteration) {
				if (bias)
					pNode = sTree.get(g1.nextInt(java.lang.Math.min(sTree.size(), maxSampleSize)));
				else
					pNode = sTree.get(g1.nextInt(sTree.size()));
				mol = newRandomMolecule(pNode.getMolecule());
				sv = mol.getRMSD(target);
				if (nearestDistance(sTree, sv)<closeDistance) {
					sIter++;
					continue;
				}
				if (rTree.size()>0 && nearestDistance(rTree, sv)<closeDistance) {
					rIter++;
					continue;
				}
				break;
			}
			if (rIter>=maxFailIteration || sIter>=maxFailIteration) {
				//System.out.println("New random configuration is too close to already sampled ones.");
				failNum++;
				continue;
			}
			
			double score = sv;
			// accept this as new branch direction based on MC criterion
			//double oldScore = minScore;
			double oldScore = pNode.score;
			if (/*bias && */score >= oldScore &&
					g2.nextDouble() >= java.lang.Math.exp((oldScore-score)/(score*acceptConstant))) {
				//System.out.println("New random configuration is not getting close to goal.");
				failNum++;
				continue;
			}
			
			// new random molecule is now accepted to make a new branch towards it
			rTree.add(new RRTNode(null, sv, null, 0));
			RRTNode newNode = extendRRT(treePool, new RRTNode(mol, sv), false);
			if (newNode == null) {
				failNum++;
				continue;
			}
			
			double newScore = newNode.score;
			int i = 0;
			for (RRTNode n:sTree) {
				if (n.score < newScore) i++;
				if (i >= maxSampleSize) break;
			}
			if (i < maxSampleSize) sTree.add(i, newNode);
			else sTree.add(newNode);
			
			failNum = 0;iter++;
			//System.out.println("Iteration: "+iter+", New Score: "+newScore);
			if (iter%20 == 0) System.out.printf(" %d", iter);
			if (iter%2000 == 0) System.out.println();
			
			BBMolecule newMol = newNode.getMolecule();
			logger.info(String.format("%.2f", newScore));
			DCDIO.writeData(out, newMol);
			if (newScore < scoreThreshold) {
				sLinkNode = newNode;tLinkNode = tNode;
				DCDIO.closeOutHandle(out);
				DCDIO.writeFrameNumber(base+dcdFile, (iter+1));
				return getPath();
			} else if (newScore < minScore) minScore = newScore;
			//bias = !bias;	// bias every other step
		}
		DCDIO.closeOutHandle(out);
		DCDIO.writeFrameNumber(base+dcdFile, iter);
		return null;
	}
	
	/**
	 * Get the path from source to target conformation 
	 */
	private ArrayList<BBMolecule> getPath() {
		if (sLinkNode == null || tLinkNode == null) return null;

		ArrayList<RRTNode> sPath = getPath(sLinkNode, sNode);
		ArrayList<RRTNode> tPath = getPath(tLinkNode, tNode);
		if (sPath == null || tPath == null) return null;
		
		log += "Source tree size: "+sTree.size()+"\n";
		log += "Target tree size: "+tTree.size()+"\n";
		//int num = sNode.getMolecule().getElementList().size();
		log += String.format(" Distance threshold: %.2f\n",  scoreThreshold);
		log += "\n";
		ArrayList<BBMolecule> path = new ArrayList<BBMolecule>();
		int size = sPath.size();
		path.add(sPath.get(size-1).getMolecule());	//add source molecule
		for (int i=size-2;i>=0;i--) {
			RRTNode n = sPath.get(i);
			path.add(n.getMolecule());
		}
		double d = sLinkNode.score- tLinkNode.score;
		log += String.format("Closest distance found between source and target tree: %.2f\n", d);
		size = tPath.size();
		for (int i=0;i<size-1;i++) {
			RRTNode n = tPath.get(i);
			path.add(n.getMolecule());
		}
		path.add(tPath.get(size-1).getMolecule());	//add target molecule
		return path;
	}

	private double getSquareDistance(RRTNode n1, RRTNode n2) {
		if (n1 == null || n2 == null) return -1;
		double sv1 = n1.score, sv2 = n2.score;
		
                
		return Math.abs(sv1-sv2);
	}
//check
	public ArrayList<RRTNode> getPath(RRTNode source, RRTNode target) {
		if (source == null || target == null) return null;
		ArrayList<RRTNode> path = null;
		if (source == target) {
			path = new ArrayList<RRTNode>();
			path.add(target);
			return path;
		}
		path = getPath(source.parent, target);
		if (path != null) {
			path.add(0, source);
			return path;
		}
		path = getPath(source, target.parent);
		if (path != null) {
			path.add(target);
			return path;
		}
		return null;
	}
        
        //check the energy if new molecule
      
	

	public String getLog() {
		return log;
	}
	
	private class RRTNode {
		private BBMolecule molecule=null;
		RRTNode parent=null;
		ArrayList<Double> step_vector=null;
		int steps = 0;
		double score ;

		/**
		 * this constructor is used by source, target and random molecules that have no parents
		 * @param mol
		 * @param scoreVector
		 */
		public RRTNode(BBMolecule mol, double score_vector) {
			molecule = mol;
			this.score = score_vector;
		}
		
		/**
		 * this constructor do not cache reference to molecule so memory can be freed
		 */
		public RRTNode(RRTNode parent,double score_vector, ArrayList<Double> step_vector, int steps) {
			this.parent = parent;
			this.score = score_vector;
			this.step_vector=step_vector;
			this.steps = steps;
		}

		/**
		 * create the molecule upon needed
		 * @return
		 */
			public BBMolecule getMolecule() {
				if (molecule != null) return molecule;
				if (parent == null) return null;
				BBMolecule mol = parent.getMolecule();
				if (mol == null) return null;
				if (step_vector == null) return null;
                                molecule = BBMolecule.newBBMolecule(mol);
				ArrayList<Integer>  l2 = mol.getDihedralAngleList();
				int  s2 = l2.size();
				for (int i=0;i<s2;i++) {
					double d = step_vector.get(i)*steps;
					// repeat what is done in extendRRT, otherwise error happens
					if (d > 0.2) {
                                            
						 molecule.rotateAroundBond(l2.get(i), l2.get(i)+1, d);
					}
				}
				RRT.this.correctAtomDistance(molecule);
				return molecule;
		}
		
	}
	
}
