using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.DrugDeliveryModel.Tests.Materials;
using MGroup.FEM.Structural.Continuum;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.Solvers.Direct;
using MGroup.MSolve.Numerics.Integration;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion;

namespace MGroup.DrugDeliveryModel.Tests.EquationModels
{
	public class Eq9ModelProviderForStaggeredSolutionex7ref
    {
        //private double sc = 0.1;
        private double miNormal;// = 5;//KPa
        private double kappaNormal;// = 6.667; //Kpa
        private double miTumor;// = 22.44; //Kpa
        private double kappaTumor;// = 216.7; //Kpa


        private List<(int, StructuralDof[], double[][], double[])> eq9LoadsList;
        private StructuralDof loadedDof; //TODO Orestis :if AddLoads9BCs() is implemented in a right way thhese will not be necessary and be deleted.
        public double load_value { get; private set; }//TODO Orestis :if AddLoads9BCs() is implemented in a right way thhese will not be necessary and be deleted.

        private List<(int, StructuralDof[], double[][], double[])> eq9BCsList;
        private double modelMinX; //TODO Orestis :if AddEquation9BCs() is implemented in a right way thhese will not be necessary and be deleted.
        private double modelMaxX;
        private double modelMinY;
        private double modelMaxY;
        private double modelMinZ;
        private double modelMaxZ;
    
        private ComsolMeshReader reader;
       
        private Dictionary<int, double> lambda;
        Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;

        public int nodeIdToMonitor { get; private set; } //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        private StructuralDof dofTypeToMonitor;
        
        
        Dictionary<int, double[]> elementslastConvergedDisplacements;
        private bool elementSavedDisplacementsIsInitialized = false;

        public int loadedNode_Id { get; private set; }

        public Eq9ModelProviderForStaggeredSolutionex7ref(ComsolMeshReader comsolReader, double sc, double miNormal, double kappaNormal, double miTumor,
            double kappaTumor, double timeStep, double totalTime,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            int nodeIdToMonitor, StructuralDof dofTypeToMonitor, StructuralDof loadedDof,
            double load_value, double modelMinX, double modelMaxX, double modelMinY, double modelMaxY, double modelMinZ, double modelMaxZ,
            List<(int, StructuralDof[], double[][], double[])> eq9BCsList, List<(int, StructuralDof[], double[][], double[])> eq9LoadsList)
        {
            //this.sc = sc;
            this.miNormal = miNormal;
            this.kappaNormal = kappaNormal;
            this.miTumor = miTumor;
            this.kappaTumor = kappaTumor;
            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            reader = comsolReader;

            this.modelMinX = modelMinX;
            this.modelMaxX = modelMaxX;
            this.modelMinY = modelMinY;
            this.modelMaxY = modelMaxY;
            this.modelMinZ = modelMinZ;
            this.modelMaxZ = modelMaxZ;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;

            //Load
            this.loadedDof = loadedDof;
            this.load_value = load_value;

            this.eq9BCsList = eq9BCsList;
            this.eq9LoadsList = eq9LoadsList;

        }

        public Model GetModel()
        {
            var nodes = reader.NodesDictionary;
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in nodes.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }
            
            var materialNormal = new NeoHookeanMaterial3dJ3Isochoric(miNormal, kappaNormal);
            var materialTumor = new NeoHookeanMaterial3dJ3Isochoric(miTumor, kappaTumor);

            var elasticMaterial = new ElasticMaterial3D(youngModulus: 1, poissonRatio: 0.3);
            var DynamicMaterial = new TransientAnalysisProperties(density: 1, rayleighCoeffMass: 0, rayleighCoeffStiffness: 0);
            var elementFactory = new ContinuumElement3DFactory(elasticMaterial, DynamicMaterial);

            //var domains = new Dictionary<int, double[]>(2);
            foreach (var elementConnectivity in reader.ElementConnectivity)
            {
                var domainId = elementConnectivity.Value.Item3;
                var element = elementFactory.CreateNonLinearElementGrowt(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2, domainId == 0 ? materialTumor : materialNormal, DynamicMaterial, lambda[elementConnectivity.Key]);
                element.volumeForce = pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0];
                element.ID = elementConnectivity.Key;
                if (elementSavedDisplacementsIsInitialized) { element.lastConvergedDisplacements = elementslastConvergedDisplacements[element.ID]; }
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }
            return model;
        }

        public void AddEquation9Loads(Model model)
        {
            foreach (var Loaddata in eq9LoadsList)
            {
                var RegionType = Loaddata.Item1; //Todo use enum suchh Region.SingleEntityPlusXFace or Region.InnerEntityMinusXFace
                var LoadedDofs = Loaddata.Item2;
                double[][] CaracteristicCoords = Loaddata.Item3;
                double[] loadsOfRespectiveDofs = Loaddata.Item4;

                switch (RegionType)
                {
                    case 1:
                        {
                            //TODO Orestis: Prosarmose tin kaloumeni parakatw methodo
                            // na min kanei xrisi kanenos field kai na 
                            // axiopoiei to caracteristic coords kai to Loaddata
                            AddEq9ModelLoadsCorner(model);
                            break;
                        }

                    case 2:
                        {
                            //TODO Orestis: (einai to validate velocity paradeigma)
                            //
                            //Omoiws
                            //  prosarmose tin 
                            // na min kanei xrisi kanenos field kai na 
                            // axiopoiei to caracteristic coords

                            //AddEq9ModelLoads(model);
                            break;
                        }

                    case 3:
                        {
                            //TODO Orestis:  Se aftin peta not Implemented exception 
                            // alla na uparxei sa methodos me swsta orismata

                            //Add_Distributed_Load(model) me oloklirwma enos double[]
                            //se mia epilegomeni epifaneia
                            break;
                        }

                }
            }

        }

        public void AddEq9ModelLoads(Model model)
        {
            INode maxDistanceNode = null;
            double currentMaxDistance = 0;
            foreach (INode node in model.NodesDictionary.Values)
            {
                double distance = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                if (distance > currentMaxDistance)
                {
                    currentMaxDistance = distance;
                    maxDistanceNode = node;
                }
            }


            /*loadedNode_Id = maxDistanceNode.ID;
            nodeIdToMonitor = loadedNode_Id;*/
            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                maxDistanceNode,
                loadedDof,
                amount: load_value
            ));

            var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));
        }

        public void AddEq9ModelLoadsCenter(Model model)
        {
            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                model.NodesDictionary[nodeIdToMonitor],
                loadedDof,
                amount: load_value
            ));

            var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));

            
        }
        
        public void AddEq9ModelLoadsCorner(Model model)
        {
            var loads = new List<INodalLoadBoundaryCondition>();

            var cornerNodes = new List<INode>();
            foreach (INode node in model.NodesDictionary.Values)
            {
                if (node.X == modelMinX && node.Y == modelMinY && node.Z == modelMaxZ)
                {
                    cornerNodes.Add(node);
                }
                if (node.X == modelMaxX && node.Y == modelMinY && node.Z == modelMaxZ)
                {
                    cornerNodes.Add(node);
                }
                if (node.X == modelMinX && node.Y == modelMaxY && node.Z == modelMaxZ)
                {
                    cornerNodes.Add(node);
                }
                if (node.X == modelMaxX && node.Y == modelMaxY && node.Z == modelMaxZ)
                {
                    cornerNodes.Add(node);
                }
            }

            foreach (INode node in cornerNodes)
            {
                loads.Add(new NodalLoad
                (
                    node,
                    loadedDof,
                    amount: load_value
                ));
            }

            var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));
        }

        public void AddEquation9BCs(Model model)
        {
            foreach (var BCdata in eq9BCsList)
            {
                var RegionType = BCdata.Item1;
                var Bcstype = BCdata.Item2;
                double[][] CaracteristicCoords = BCdata.Item3;
                double[] prescrVal = BCdata.Item4;

                switch (RegionType)
                {
                    case 1:
                        {
                            //TODO Orestis: Prosarmose tin kaloumeni parakatw methodo
                            // na min kanei xrisi kanenos field kai na 
                            // axiopoiei to caracteristic coords kai to BCdata
                            AddBottomLeftRightFrontBackBCs(model);
                            break;
                        }

                    case 2:
                        {
                            //TODO Orestis: copy  AddBottomBCs() method from
                            //Eq9ModelProviderForStaggeredSolutionex83.cs
                            // kai prosarmose tin 
                            // na min kanei xrisi kanenos field kai na 
                            // axiopoiei to caracteristic coords

                            //AddAllBoundaryNodesBC(model);
                            break;
                        }

                    case 3:
                        {
                            //TODO Orestis: copy  AddBottomLeftFrontBackBCs() method from
                            //Eq9ModelProviderForStaggeredSolutionex6_1.cs
                            // kai prosarmose tin 
                            // na min kanei xrisi kanenos field kai na 
                            // axiopoiei to caracteristic coords

                            //AddAllBoundaryNodesBC(model);
                            break;
                        }

                }
            }

        }

        public void AddXFaceBcs(Model model, double xCoordOfFace, double[] dirichletValuesToPrescribe, StructuralDof[] dofTypesToConstrained)
        {

            var xFaceNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(xCoordOfFace - node.X) < 1E-9) xFaceNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();

            foreach (var node in xFaceNodes)
            {
                for (int i = 0; i < dirichletValuesToPrescribe.Length; i++)
                {
                    var dofToSetPrescribed = dofTypesToConstrained[i];
                    var valueToPrescribe = dirichletValuesToPrescribe[i];
                    constraints.Add(new NodalDisplacement(node, dofToSetPrescribed, amount: valueToPrescribe));
                }
                
            }

            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));

        }

        public void AddYFaceBcs(Model model, double yCoordOfFace, double[] dirichletValuesToPrescribe, StructuralDof[] dofTypesToConstrained)
        {

            var yFaceNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(yCoordOfFace - node.Y) < 1E-9) yFaceNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();

            foreach (var node in yFaceNodes)
            {
                for (int i = 0; i < dirichletValuesToPrescribe.Length; i++)
                {
                    var dofToSetPrescribed = dofTypesToConstrained[i];
                    var valueToPrescribe = dirichletValuesToPrescribe[i];
                    constraints.Add(new NodalDisplacement(node, dofToSetPrescribed, amount: valueToPrescribe));
                }

            }

            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));

        }


        public void AddZFaceBcs(Model model, double zCoordOfFace, double[] dirichletValuesToPrescribe, StructuralDof[] dofTypesToConstrained)
        {

            var zFaceNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(zCoordOfFace - node.Z) < 1E-9) zFaceNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();

            foreach (var node in zFaceNodes)
            {
                for (int i = 0; i < dirichletValuesToPrescribe.Length; i++)
                {
                    var dofToSetPrescribed = dofTypesToConstrained[i];
                    var valueToPrescribe = dirichletValuesToPrescribe[i];
                    constraints.Add(new NodalDisplacement(node, dofToSetPrescribed, amount: valueToPrescribe));
                }

            }

            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));

        }
        public void AddXplusBCs(Model model, double modelMaxX, double ValueprescribedformaxXcoord)
        {

            var maxXCoordNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxX - node.X) < 1E-9) maxXCoordNodes.Add(node);
            }

            var dirichletBCs = new List<NodalUnknownVariable>();
            foreach (var node in maxXCoordNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, ValueprescribedformaxXcoord));
            }


            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

        }

        public void AddYplusBCs(Model model, double modelMaxX, double ValueprescribedformaxXcoord)
        {

            var maxXCoordNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxX - node.X) < 1E-9) maxXCoordNodes.Add(node);
            }

            var dirichletBCs = new List<NodalUnknownVariable>();
            foreach (var node in maxXCoordNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, ValueprescribedformaxXcoord));
            }


            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

        }

        public void AddBottomLeftRightFrontBackBCs(Model model)
        {
            var bottomNodes = new List<INode>();
            var leftNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var frontNodes = new List<INode>();
            var backNodes = new List<INode>();

            var innerBulkNodes = new List<INode>();
            var tol = 1E-5;
            //Search for all boundary nodes
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) bottomNodes.Add(node);
                if(Math.Abs(modelMinX - node.X) < tol) leftNodes.Add(node);
                if(Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);
                if(Math.Abs(modelMinY - node.Y) < tol) frontNodes.Add(node);
                if(Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
            }

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) { }

                else if (Math.Abs(modelMinX - node.X) < tol) { } 
                else if (Math.Abs(modelMaxX - node.X) < tol) { }

                else if (Math.Abs(modelMinY - node.Y) < tol) { }
                else if (Math.Abs(modelMaxY - node.Y) < tol) { }
                else innerBulkNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();
            
            //Apply roller constraint to bottom nodes. (constrained movement in z direction)
            foreach (var node in bottomNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));

            //Apply roller constraint to left nodes. (constrained movement in x direction)
            foreach (var node in leftNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            
            //Apply roller constraint to right nodes. (constrained movement in x direction)
            foreach (var node in rightNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            
            //Apply roller constraint to front nodes. (constrained movement in y direction)
            foreach (var node in frontNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            
            //Apply roller constraint to back nodes. (constrained movement in y direction)
            foreach (var node in backNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            
            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));
        }
        
        //TODO Gerasimos add if for dynamic or peudostatic analyzer

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel =  solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemStructural(model, algebraicModel);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: nIncrements)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var nlAnalyzer = loadControlAnalyzerBuilder.Build();
            var loadControlAnalyzer = (LoadControlAnalyzer)nlAnalyzer;
            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    

                }, algebraicModel
            );

            //var analyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: currentStep)).Build();
            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime,false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            //var analyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentTimeStep: currentStep, bdfOrder: 5);
            var analyzer = analyzerBuilder.Build();


            //Sparse tet Mesh
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    //(model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationY),
                    //(model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationZ),
                }
            };

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }

        public void SaveStateFromElements(Model model)
        {
            elementslastConvergedDisplacements = new Dictionary<int, double[]>();
            foreach (var elem in reader.ElementConnectivity)
            {
                elementslastConvergedDisplacements[elem.Key] = ((ContinuumElement3DGrowth)model.ElementsDictionary[elem.Key]).localDisplacements.Copy();
            }

            elementSavedDisplacementsIsInitialized = true;
        }
    }
}
