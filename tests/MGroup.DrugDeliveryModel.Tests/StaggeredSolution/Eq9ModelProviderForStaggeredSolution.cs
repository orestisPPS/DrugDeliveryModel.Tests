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

namespace MGroup.DrugDeliveryModel.Tests.EquationModels
{
	public class Eq9ModelProviderForStaggeredSolution
    {
       

        //private double sc = 0.1;
        private double miNormal;// = 5;//KPa
        private double kappaNormal;// = 6.667; //Kpa
        private double miTumor;// = 22.44; //Kpa
        private double kappaTumor;// = 216.7; //Kpa
        private StructuralDof loadedDof;
        private double load_value;


        private ComsolMeshReader reader;
       
        private Dictionary<int, double> lambda;
        Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;

        private int nodeIdToMonitor; //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        private StructuralDof dofTypeToMonitor;
        

        private double modelMaxZ;
        private double topValueprescribed;
        private double modelMinZ;
        private double bottomValueprescribed;




        public Eq9ModelProviderForStaggeredSolution(ComsolMeshReader comsolReader, double sc, double miNormal, double kappaNormal, double miTumor,
            double kappaTumor, double timeStep, double totalTime,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            double modelMaxZ, double topValueprescribed, double modelMinZ, double bottomValueprescribed,
            int nodeIdToMonitor, StructuralDof dofTypeToMonitor, StructuralDof loadedDof, double load_value)
        {
            //this.sc = sc;
            this.miNormal = miNormal;
            this.kappaNormal = kappaNormal;
            this.miTumor = miTumor;
            this.kappaTumor = kappaTumor;
            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            reader = comsolReader;

            //BCs 1
            this.modelMaxZ = modelMaxZ;
            this.topValueprescribed = topValueprescribed;
            this.modelMinZ = modelMinZ;
            this.bottomValueprescribed = bottomValueprescribed;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;

            //Load
            this.loadedDof = loadedDof;
            this.load_value = load_value;

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

            //(double ENormal, double poissonNormal) = GetElasticModelParameters(miNormal, kappaNormal);
            //(double ETumor, double poissonTumor) = GetElasticModelParameters(miTumor, kappaTumor);
            //var materialNormal = new ElasticMaterial3DDefGrad(ENormal, poissonNormal);
            //var materialTumor = new ElasticMaterial3DDefGrad(ETumor, poissonTumor);
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
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

            return model;
        }

        public void AddEq9ModelAppropriateBCs(Model model)
        {
            var faceXYNodes = new List<INode>();
            var faceXZNodes = new List<INode>();
            var faceYZNodes = new List<INode>();

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(0 - node.Z) < 1E-9) faceXYNodes.Add(node);
                if (Math.Abs(0 - node.Y) < 1E-9) faceXZNodes.Add(node);
                if (Math.Abs(0 - node.X) < 1E-9) faceYZNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();
            foreach (var node in faceXYNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));
            }
            foreach (var node in faceXZNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            }
            foreach (var node in faceYZNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            }

            var emptyLoads1 = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyLoads1));

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

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel =  solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemStructural(model, algebraicModel);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: nIncrements)
            {
                ResidualTolerance = 1E-4,
                MaxIterationsPerIncrement = 1000,
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

            var analyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: currentStep)).Build();

            //Sparse tet Mesh
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    //(model[0].NodesDictionary[333], StructuralDof.TranslationY),
                    //(model[0].NodesDictionary[333], StructuralDof.TranslationZ),
                }
            };


            //Print the coordinates of the watchdof for comsol comparison
            //Console.WriteLine("DOF :" + model[0].NodesDictionary[13].ID + " X: " + model[0].NodesDictionary[13].X + " Y: " + model[0].NodesDictionary[13].Y + " Z: " + model[0].NodesDictionary[13].Z + " Time: " + currentTimeStep * timeStep + " lambda: " + lambda[13]);
            //Console.WriteLine("DOF :" + model[0].NodesDictionary[333].ID + " X: " + model[0].NodesDictionary[333].X + " Y: " + model[0].NodesDictionary[333].Y + " Z: " + model[0].NodesDictionary[333].Z + " Time: " + currentTimeStep * timeStep + " lambda: " + lambda[333]);

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }

        //public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        //{
        //    //Changed this
            

        //    var model = new Model[] { CreateElasticModelFromComsolFile(reader, miNormal, kappaNormal, miTumor, kappaTumor, lambda), };
        //    var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
        //    var algebraicModel = new[] { solverFactory.BuildAlgebraicModel(model[0]), };
        //    solvers[0] = solverFactory.BuildSolver(algebraicModel[0]);
        //    var problem = new[] { new ProblemStructural(model[0], algebraicModel[0]), };
        //    var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder( algebraicModel[0], solvers[0], problem[0], numIncrements: nrIncrements)
        //    {
        //        ResidualTolerance = 1E-4,
        //        MaxIterationsPerIncrement = 1000,
        //        NumIterationsForMatrixRebuild = 1
        //    };
        //    nlAnalyzers[0] = loadControlAnalyzerBuilder.Build();
        //    var loadControlAnalyzer = (LoadControlAnalyzer)nlAnalyzers[0];
        //    loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
        //        new List<(INode node, IDofType dof)>()
        //        {
        //            (model[0].NodesDictionary[333], StructuralDof.TranslationX),
        //            (model[0].NodesDictionary[333], StructuralDof.TranslationY),
        //            (model[0].NodesDictionary[333], StructuralDof.TranslationZ),

        //        }, algebraicModel[0]
        //    );

        //    analyzers[0] = (new PseudoTransientAnalyzer.Builder( algebraicModel[0], problem[0], loadControlAnalyzer, timeStep: timeStep, totalTime: totalTime, currentStep: CurrentTimeStep)).Build();

        //    //Sparse tet Mesh
        //    var watchDofs = new[]
        //    {
        //        new List<(INode node, IDofType dof)>()
        //        {
        //            (model[0].NodesDictionary[333], StructuralDof.TranslationX),
        //            (model[0].NodesDictionary[333], StructuralDof.TranslationY),
        //            (model[0].NodesDictionary[333], StructuralDof.TranslationZ),
        //        }
        //    };


        //    //Print the coordinates of the watchdof for comsol comparison
        //    //Console.WriteLine("DOF :" + model[0].NodesDictionary[13].ID + " X: " + model[0].NodesDictionary[13].X + " Y: " + model[0].NodesDictionary[13].Y + " Z: " + model[0].NodesDictionary[13].Z + " Time: " + currentTimeStep * timeStep + " lambda: " + lambda[13]);
        //    //Console.WriteLine("DOF :" + model[0].NodesDictionary[333].ID + " X: " + model[0].NodesDictionary[333].X + " Y: " + model[0].NodesDictionary[333].Y + " Z: " + model[0].NodesDictionary[333].Z + " Time: " + currentTimeStep * timeStep + " lambda: " + lambda[333]);

        //    loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel[0]);

        //    for (int i = 0; i < analyzers.Length; i++)
        //    {
        //        analyzers[i].Initialize(true);
        //        if (analyzerStates[i] != null)
        //        {
        //            analyzers[i].CurrentState = analyzerStates[i];
        //        }

        //        if (nlAnalyzerStates[i] != null)
        //        {
        //            nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
        //        }
        //    }
        //}

        public static Model CreateElasticModelFromComsolFile(ComsolMeshReader reader, double miNormal, double kappaNormal, double miTumor, double kappaTumor, Dictionary<int, double> lambda)
        {
            var nodes = reader.NodesDictionary;
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in nodes.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }

            //(double ENormal, double poissonNormal) = GetElasticModelParameters(miNormal, kappaNormal);
            //(double ETumor, double poissonTumor) = GetElasticModelParameters(miTumor, kappaTumor);
            //var materialNormal = new ElasticMaterial3DDefGrad(ENormal, poissonNormal);
            //var materialTumor = new ElasticMaterial3DDefGrad(ETumor, poissonTumor);
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
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

            var faceXYNodes = new List<INode>();
            var faceXZNodes = new List<INode>();
            var faceYZNodes = new List<INode>();

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(0 - node.Z) < 1E-9) faceXYNodes.Add(node);
                if (Math.Abs(0 - node.Y) < 1E-9) faceXZNodes.Add(node);
                if (Math.Abs(0 - node.X) < 1E-9) faceYZNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();
            foreach (var node in faceXYNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));
            }
            foreach (var node in faceXZNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
            }
            foreach (var node in faceYZNodes)
            {
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
            }

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


            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                maxDistanceNode,
                StructuralDof.TranslationX,
                amount: 0.00000001
            ));


            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

            return model;
        }

        /// <summary>
        /// https://en.wikiversity.org/wiki/Elasticity/Constitutive_relations
        /// </summary>
        /// <param name="mi"></param>
        /// <param name="kappa"></param>
        /// <returns></returns>
        static (double, double) GetElasticModelParameters(double mi, double kappa)
        {
            double poisson = (3 * kappa - 2 * mi) / (2 * (3 * kappa + mi));
            double E = 9 * kappa * mi / (3 * kappa + mi);
            return (E, poisson);
        }
    }
}
