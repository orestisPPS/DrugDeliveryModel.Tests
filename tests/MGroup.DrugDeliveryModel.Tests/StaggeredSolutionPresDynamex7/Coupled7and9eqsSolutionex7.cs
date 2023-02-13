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
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled7and9eqsSolutionex7
    {
        const double Sc = 0.1;

        private const double timeStep = 0.00001; // in sec
        const double totalTime = 0.001; // in sec
        static int incrementsPertimeStep = 1;


        // strucutral model Loads
        static StructuralDof loadedDof = StructuralDof.TranslationZ;
        //static double load_value = 0.0001 / 4; //[kN]
        static double load_value = 0.0001 / 4; //[kN]


        //structural model properties
        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa

        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa

        static int currentTimeStep = 0;
        static double lambda0 = 1;
        int nGaussPoints = 1;
        static double initial_dp_dx = 0.0;
        static double initial_dp_dy = 0.0;
        static double initial_dp_dz = 0.0;
        static double velocityDivInitialVal = 0;
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>();
        private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();


        private static int pressureMonitorID;

        /*static double divPMonitorGPX = 0.01714618146;
        static double divPMonitorGPY = 0.05114110504;
        static double divPMonitorGPZ = 0.06533804930;*/

        //static double divPMonitorGPX = 0.05;
        //static double divPMonitorGPY = 0.05;
        //static double divPMonitorGPZ = 0.05;

        //static double divPMonitorGPX = 0;
        //static double divPMonitorGPY = 0;
        //static double divPMonitorGPZ = 0.1;

        static double divPMonitorGPX = 0.004086132769345323;
        static double divPMonitorGPY = 0.006191771651571988;
        static double divPMonitorGPZ = 0.09369050147682026;
                                       
        //GP with coordinates: 0.05301208792514899 0.053825572057669926 0.052065045951539365 is in element with id: 1251 NEAR [0.5, 0.5, 0.5]
        //GP with coordinates: 0.004086132769345323 0.006191771651571988 0.09369050147682026 is in element with id: 740  NEAR [0,0,0.1]


        //static double[] monitoredGPcoords = new double[] { 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 };
        static double[] monitoredGPcoords = new double[] { 0.004086132769345323, 0.006191771651571988, 0.09369050147682026 };

        private static int divPMonitorID;

        static StructuralDof eq9dofTypeToMonitor = StructuralDof.TranslationZ;

        //Darcy model properties
        
        // original case
        //static double Sv = 7e+3; // 1/(m)
        //static double Lp = 2.7e-12; // m/(KPa sec)
        //static double LplSvl = 0; // 1/(KPa sec)
        //static double pv = 4; // kPa
        //static double pl = 0d; // KPa
        //static double k_th = 7.52e-10; // m2/(KPa sec)

        // simplified case
        static double Sv = 0; // 1/(m)
        static double Lp = 0; // m/(KPa sec)
        static double LplSvl = 0; // 1/(KPa sec)
        static double pv = 0; // kPa
        static double pl = 0d; // KPa
        static double k_th = 7.52e-6; // m2/(KPa sec)




        //static int nodeIdToMonitor = 36;

        //The above coordinates are in the reference configuration.
        //The actual coordinates will be computed in the analysis.

        //static double structuralMonitorNodeX = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa
        //static double structuralMonitorNodeY = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa
        //static double structuralMonitorNodeZ = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa

        static double structuralMonitorNodeX = 0.0;
        static double structuralMonitorNodeY = 0.0;
        static double structuralMonitorNodeZ = 0.1;

        //static double[] structuralMonitorNodeCoords = new double[]
        //    { 0.04930793848882013,0.04994681648346263,0.04953188199244812 };
        //static double[] structuralMonitorNodeCoords = new double[]
        //    { 0.04379781548934232, 0.0492383785341323, 0.02519415945383001 };
        static double[] structuralMonitorNodeCoords = new double[]
            { 0.0, 0.0, 0.1 };

        private static int structuralMonitorID;
        static ConvectionDiffusionDof eq7n8dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;

        //private static double boundaryValueAllBoundaries = pv;
        //private static double initialCondition = pv;
        private static double boundaryValueAllBoundaries = 0;
        private static double initialCondition = 0;

        static double modelMinX = 0;
        static double modelMaxX = 0.1;
        static double modelMinY = 0;
        static double modelMaxY = 0.1;
        static double modelMinZ = 0;
        static double modelMaxZ = 0.1;

        public Coupled7and9eqsSolutionex7()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]
        //[InlineData("../../../DataFiles/workingTetMesh4886.mphtxt")]
        //[InlineData("../../../DataFiles/chipMelter2M.mphtxt")]
        //[InlineData("../../../DataFiles/MeshCyprusTM.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh4886_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh648_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh648_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingQuadMesh64_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingQHexaMesh6x6x6_1Domain.mphtxt")]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingQuadMesh27_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh155.mphtxt")]
        public void MonophasicEquationModel(string fileName)
        {
            ContinuumElement3DGrowth.dT = timeStep;

            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);

            // initialize Shared quantities of Coupled model
            Dictionary<int, double> lambda = new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                lambda.Add(elem.Key, lambda0);
            }

            Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints =
                new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var gpTensorDiv = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    gpTensorDiv[i1] = new double[] { initial_dp_dx, initial_dp_dy, initial_dp_dz };
                }

                pressureTensorDivergenceAtElementGaussPoints.Add(elem.Key, gpTensorDiv);
            }

            Dictionary<int, double[]> velocityDivergenceAtElementGaussPoints =
                new Dictionary<int, double[]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocityDiv = new double[nGaussPoints];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    velocityDiv[i1] = velocityDivInitialVal;
                }

                velocityDivergenceAtElementGaussPoints.Add(elem.Key, velocityDiv);
            }

            miNormal = miTumor;
            kappaNormal = kappaTumor;

            #region loggin (defined before model builder creation to give them nodes)
            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNodeCoords, 1e-2);
            //pressureMonitorID = Utilities.FindRandomInternalNode(comsolReader.NodesDictionary, modelMinX, modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            pressureMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, new double[] { 0.055, 0.0559, 0.07366 }, 1e-2);
            //pressureMonitorID = structuralMonitorID;

            var p_i = new double[(int)(totalTime / timeStep)];

            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>();
            displacements.Add(structuralResultsX);
            displacements.Add(structuralResultsY);
            displacements.Add(structuralResultsZ);

            double[] gp_dut_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dvt_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dwt_dz_OverTime = new double[(int)(totalTime / timeStep)];
            var divVelocity = new List<double[]>();
            divVelocity.Add(gp_dut_dx_OverTime);
            divVelocity.Add(gp_dvt_dy_OverTime);
            divVelocity.Add(gp_dwt_dz_OverTime);


            double[] gp_div_v_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dz_Overtime = new double[(int)(totalTime / timeStep)];
            var dp_dxi = new List<double[]>();
            dp_dxi.Add(gp_dP_dx_OverTime);
            dp_dxi.Add(gp_dP_dy_OverTime);
            dp_dxi.Add(gp_dP_dz_Overtime);

            int monitoredGP_elemID = -1;
            #endregion




            var eq78Model = new Eq78ModelProviderForStaggeredSolutionex7(comsolReader, k_th, Lp, Sv, pv, LplSvl, pl,
                velocityDivergenceAtElementGaussPoints,

                boundaryValueAllBoundaries, initialCondition, pressureMonitorID, eq7n8dofTypeToMonitor, modelMinX,
                modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            //COMMITED BY NACHO 
            //jkkk bn///////vji typ[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[00u-----------------------------------
            var eq9Model = new Eq9ModelProviderForStaggeredSolutionex7(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, eq9dofTypeToMonitor, loadedDof, load_value, modelMinX, modelMaxX, modelMinY,
                modelMaxY, modelMinZ, modelMaxZ);

            var equationModel = new Coupled7and9eqsModelex7(eq78Model, eq9Model, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, timeStep,
                totalTime, incrementsPertimeStep);



            


            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.001);                                                       
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging
                monitoredGP_elemID =Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoords, 1e-1);

                //nodal logs
                p_i[currentTimeStep] =((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), ConvectionDiffusionDof.UnknownVariable];
                //p_i[currentTimeStep] = 0d;
                //structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                //structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                //structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];
                structuralResultsX[currentTimeStep] = 0d;
                structuralResultsY[currentTimeStep] = 0d;
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];

                //gp (element) logs
                gp_dP_dx_OverTime[currentTimeStep] =((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGP_elemID]).xcoeff_OverTimeAtGp1[0];
                gp_dP_dy_OverTime[currentTimeStep] =((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGP_elemID]).ycoeff_OverTimeAtGp1[0];
                gp_dP_dz_Overtime[currentTimeStep] =((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGP_elemID]).zcoeff_OverTimeAtGp1[0];
                gp_dut_dx_OverTime[currentTimeStep]= ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence_term1[0];
                gp_dvt_dy_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence_term2[0];
                gp_dwt_dz_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence_term3[0];
                gp_div_v_OverTime[currentTimeStep]= ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence[0];

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

                (equationModel.ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
                (equationModel.ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                equationModel.SaveStateFromElements();

                //Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            //Assert.True(ResultChecker.CheckResults(structuralResultsZ, expectedDisplacments(), 1E-6));
            //Assert.True(ResultChecker.CheckResults(structuralResultsZ, expectedPressurevalues(), 1E-6));

            double pr = 100;
            double Fval = 100;
            double Fval_e = 0;
            //eq9model.load_value;

            //pressure_{ pr}_F_{Fval}_e{Fval_e}_LOGGEDval_ PAth name do not erase this exampleNo_{exNo}_caseNo_{caseNo}_
            var writer = new MGroup.LinearAlgebra.Output.Array1DWriter();
            var patSelection = 0;
            var outputPath = patSelection == 0 ? "../../../StaggeredSolutionPresDynamex7/results1/" : $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\";
            writer.WriteToFile(structuralResultsX, outputPath + $@"u_{1}_.txt");
            writer.WriteToFile(structuralResultsY, outputPath + $@"u_{2}_.txt");
            writer.WriteToFile(structuralResultsZ, outputPath + $@"u_{3}_.txt");
            writer.WriteToFile(gp_dut_dx_OverTime, outputPath + $@"gp_dut_dx_OverTime.txt");
            writer.WriteToFile(gp_dvt_dy_OverTime, outputPath + $@"gp_dvt_dy_OverTime.txt");
            writer.WriteToFile(gp_dwt_dz_OverTime, outputPath + $@"gp_dwt_dz_OverTime.txt");
            writer.WriteToFile(gp_div_v_OverTime, outputPath +  $@"gp_div_v_OverTime_.txt");



            writer.WriteToFile(p_i, outputPath + $@"p_i.txt");
            writer.WriteToFile(gp_dP_dx_OverTime, outputPath + $@"gp_dP_dx_OverTime.txt");
            writer.WriteToFile(gp_dP_dy_OverTime, outputPath + $@"gp_dP_dy_OverTime.txt");
            writer.WriteToFile(gp_dP_dz_Overtime, outputPath + $@"gp_dP_dz_Overtime.txt");


            //var path = outputPath+"dp_dxi_mslv.csv";
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../StaggeredSolutionPresDynamex7/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../StaggeredSolutionPresDynamex7/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../StaggeredSolutionPresDynamex7/pi_nodes_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(divVelocity), "../../../StaggeredSolutionPresDynamex7/dut_dxi_GP_mslv.csv");

        }


        public static double[] expectedDisplacments()
        {
            return new double[] { };
        }


        public static double[] expectedPressurevalues()
        {
            return new double[] { };
        }






    }
}
