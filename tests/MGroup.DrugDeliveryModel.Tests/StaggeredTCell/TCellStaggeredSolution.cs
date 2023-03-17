using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using TriangleNet.Meshing.Algorithm;
using System.Xml.Linq;
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using BC = MGroup.DrugDeliveryModel.Tests.Commons.BoundaryAndInitialConditionsUtility.BoundaryConditionCase;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class TCellStaggeredSolution
    {
        const double Sc = 0.1;

        private const double timeStep = 1E-5; // in sec
        const double totalTime = 1E-2; // in sec
        static int incrementsPertimeStep = 1;
        static int currentTimeStep = 0;
        
        #region Cancer Cell Density (TCell) model

        private const double dummySolidVelovity = -5;

        /// <summary>
        /// Growth rate parameter 1[mol/(m3)]
        /// </summary>
        private const double K1 = 1.74E-6; // [1/s]

        /// <summary>
        /// Growth rate parameter 2[mol/(m3)]
        /// </summary>
        private const double K2 = 8.3E-3; // [mol/(m3)]

        /// <summary>
        /// Oxygen concentration (Dependent Variable) [mol/m3]
        /// </summary>
        private const double Cox = 0.2; // [mol / m3]

        #endregion

        #region Cancer Cell Density (TCell) Boundary Conditions

        private static ConvectionDiffusionDof[] constrainedDofType = new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable };
        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])> tCellDirichletBC =
            new List<(BC, ConvectionDiffusionDof[], double[][], double[])>()
        {(BC.TopRightBackDiriclet, constrainedDofType, new double[3][]{new double[3] {0.1,0.1,0.1},new double[3] {0.1,0.1,0.1}, new double[3]{0.1,0.1,0.1}}, new double[]{500d}),};

        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])>
                tCellNeumannBC = new List<(BC, ConvectionDiffusionDof[], double[][], double[])>();

        #endregion

        #region Cancer Cell Density (TCell)  Initial condition

        private double initialTCellDensity = 0d;

        #endregion

        #region Cancer Cell Density (TCell)  logs

        static List<(double[], string, StructuralDof, int, double[])> nodeTCellLogs = new List<(double[], string, StructuralDof, int, double[])>()
            {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        //static double[] tCellMonitorNodeCoords = new double[] { 0.055, 0.0559, 0.07366 };
        static double[] tCellMonitorNodeCoords = { 0.0, 0.09, 0.09 };

        private static int tCellMonitorID;

        static ConvectionDiffusionDof tCellMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        #endregion

        public TCellStaggeredSolution()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]

        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]

        public void MonophasicEquationModel(string fileName)
        {
            ContinuumElement3DGrowth.dT = timeStep;

            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);



            Dictionary<int, double[][]> velocityAtGaussPoints =
                new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            var nGaussPoints = 1;
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocityDiv = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    velocityDiv[i1] = new double[] { 0, 0, 0 };
                };
                velocityAtGaussPoints.Add(elem.Key, velocityDiv);
            }

            Dictionary<int, double> dummyFieldCOx =
                new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                //dummyFieldCOx.Add(elem.Key, Cox);
                dummyFieldCOx.Add(elem.Key, 0d);
            }



            #region loggin (defined before model builder creation to give them nodes)
            
            tCellMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, tCellMonitorNodeCoords, 1e-2);
            
            double[] tCell = new double[(int)(totalTime / timeStep) + 1];

            #endregion
            
            //Create Model For TCell
            var tCellModel = new TCellModelProvider(K1, K2, dummyFieldCOx, velocityAtGaussPoints, comsolReader, tCellMonitorDOF, tCellMonitorID, tCellDirichletBC, tCellNeumannBC, initialTCellDensity);



            var equationModel = new TCellStaggeredModelProvider(tCellModel, comsolReader, velocityAtGaussPoints,
                timeStep, totalTime, 10);

            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.000000001);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging
                tCell[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(tCellMonitorID), tCellMonitorDOF];

                //model maximus (DO NOT ERASE)
                //modelMaxVelDivOverTime[currentTimeStep] = velocityDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0])).ToArray().Max();
                //modelMax_dP_dxOverTime[currentTimeStep] = pressureTensorDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0][0])).ToArray().Max();


                /*if (Solution.ContainsKey(currentTimeStep))
                {
                    Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    Solution.Add(currentTimeStep, allValues);
                }*/

                #endregion
                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    (equationModel.ParentAnalyzers[j] as NewmarkDynamicAnalyzer).AdvanceStep();
                }

                

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

            }
            
            //Assert.True(ResultChecker.CheckResults(tCell, expected_Tc_values(), 1e-1));
            
            CSVExporter.ExportVectorToCSV(tCell, "../../../StaggeredTCell/tCell_nodes_mslv.csv");


        }

        

        public static double[] expected_Tc_values()
        {
            return new double[] {
            0.39161779623154619,
            0.77197862637823178,
            1.14114595438075,
            1.4991842343623834,
            1.8461589027996081,
            2.1821363705976768,
            2.507184015072121,
            2.8213701718371418,
            3.1247641266018533,
            3.4174361068753587,
            };
        }







    }
}
