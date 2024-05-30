import os
import pathlib
from image_proc_functions import pheno_filt
import csv
import matplotlib 
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import random
from shapely.geometry import Polygon
from scipy.spatial.distance import cdist

def add_celltype(df, celltype_list, pheno_list):
    '''
    Function to associate cell type label to Phenotype(s)
    Args:
        df : DataFrame
        celltype_list : DataFrame
        pheno_list : list

    Returns:
        df: DataFrame
    '''
    pheno_df=pheno_filt(df, pheno_list)
    pheno_df=pheno_df[['Cell.ID', 'Cell.X.Position', 'Cell.Y.Position', 'Pheno']]
    pheno_df=pheno_df.loc[pheno_df['Pheno'].isin(celltype_list["Phenotype"])]
    
    cell_type=celltype_list.set_index('Phenotype').to_dict('dict')['Cell_Type']
    pheno_df['Celltype']=pheno_df['Pheno'].replace(cell_type)
    
    return pheno_df.reset_index()

def create_output_csv(output_path_csv, C_1, C_2):
    '''
    Function to create tsv header for summary file
    Args:
        output_path_csv : string
        C_1 : string
        C_2 : string

    Returns:
        csv_path: string
    '''
    pathlib.Path(output_path_csv).mkdir(parents=True, exist_ok=True)

    csv_path=os.path.join(output_path_csv, f'{C_1}-{C_2}.tsv')
    # Creare il file CSV e scrivere l'intestazione
    with open(csv_path, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        header = ['Subject', 'ROI', 'PCF_r', 'Counts_C1', 'Counts_C2']
        writer.writerow(header)

    return csv_path

def load_point_cloud(df, output_path):
    '''
    Function to generate and plot pointcloud
    Args:
        df : DataFrame
        output_path : string

    Returns:
        pc: pointcloud
    '''
    points = np.asarray([df['Cell.X.Position'],df['Cell.Y.Position']]).transpose()
    
    max_x=df['Cell.X.Position'].max() + 100
    max_y=df['Cell.Y.Position'].max() + 100
    min_x=df['Cell.X.Position'].min() - 100
    min_y=df['Cell.Y.Position'].min() - 100

    pc = generatePointCloud('ROI',points,domain=[[min_x,max_x],[min_y,max_y]])

    pc.addLabels('Celltype','categorical',df.Celltype,cmap='tab20')

    for i, label in enumerate(df['Celltype']):
        try:
            pc.changeIndividualLabelColor('Celltype', label, plt.cm.tab20(i))
        except:
            print(f'label {label} not found!')
            pass

    visualisePointCloud(pc, 'Celltype', cmap='tab20',markerSize=100)

    plt.savefig(os.path.join( output_path, f"point_cloud.tif"), dpi=300, format="tiff", bbox_inches='tight')
    plt.close()
    
    return pc

def all_cross_pcf(pc, output_path, roi_name, maxR=150, annulusStep=10, annulusWidth=10):
    '''
    Function to evaluate and plot all pcf combination
    Args:
        pc : pointcloud
        output_path : string
        roi_name : string
        maxR : int, optional (Defaults to 150)
        annulusStep : int, optional (Defaults to 10)
        annulusWidth : int, optional (Defaults to 10)

    Returns:
        pc: pointcloud
    '''   
    outpath=os.path.join(output_path, "All_PCF")
    pathlib.Path(outpath).mkdir(parents=True, exist_ok=True)
    #calculation of all PCFs and plotting them
    avals = []
    bvals = []
    pcfs = []
    for a in pc.labels['Celltype']['categories']:
        for b in pc.labels['Celltype']['categories']:
            #logger.info(a + " - " +b)
            r, pcf,_ = pairCorrelationFunction(pc, 'Celltype', [a,b], maxR=maxR,annulusStep=annulusStep,annulusWidth=annulusWidth)
            avals.append(a)
            bvals.append(b)
            pcfs.append(pcf)
    # the pcf threshold = 1 is caused by CSR (complete spatial randomness)
    sns.set(font_scale=1.0)
    fig, ax = plt.subplots(nrows=len(pc.labels['Celltype']['categories']), ncols=len(pc.labels['Celltype']['categories']),sharex=True,sharey=True, figsize=(20,20))
    it = 0
    for i,a in enumerate(pc.labels['Celltype']['categories']):
        for j,b in enumerate(pc.labels['Celltype']['categories']):
            ax[i,j].plot(r,pcfs[it],lw=3)
            if i == len(pc.labels['Celltype']['categories'])-1:
                ax[i,j].set_xlabel(b, rotation=45)
            if j == 0:
                ax[i,j].set_ylabel(a, rotation=45, labelpad=30)
                #ax[i, j].yaxis.set_label_coords(-0.2, 0.5)
            ax[i,j].set_ylim([0,10])
            ax[i,j].axhline(1,linestyle=':',c='k',lw=3)
            it = it + 1
    plt.savefig(os.path.join(outpath, f"ROI_{roi_name}.tif"), dpi=300, format="tiff", bbox_inches='tight')
    plt.close()

    return pc

def selected_PCF(C_1, C_2, pc, output_path, radiusOfInterest, maxR=150, annulusStep=10, annulusWidth=10):
    '''
    Function to evaluate and plot a specific pcf combination 
    Args:
        C_1 : string
        C_2 : string
        pc : pointcloud
        output_path : string
        radiusOfInterest : int
        maxR : int, optional (Defaults to 150)
        annulusStep : int, optional (Defaults to 10)
        annulusWidth : int, optional (Defaults to 10)

    Returns:
        pc : pointcloud
        pcf_value_at_radius : float

    '''
    sns.set(font_scale=2)
    
    r, pcf,_ = pairCorrelationFunction(pc, 'Celltype', [C_1,C_2], maxR=maxR,annulusStep=annulusStep,annulusWidth=annulusWidth)

    plt.figure(figsize=(15,15))
    plt.plot(r, pcf, lw=5)
    plt.gca().axhline(1,c='k',linestyle=':',lw=3)
    plt.xlabel('$r$ ($\mu$m)')
    plt.ylabel('$g_{ThE}(r)$')
    plt.figtext(0.5, 0.95, f"PCF function: {C_1} - {C_2}", ha='center', fontsize=20, bbox=dict(facecolor='white', alpha=0.5))

    plt.savefig(os.path.join(output_path, f"PCF_function.tif"), dpi=300, format="tiff", bbox_inches='tight')
    plt.close()

    # extracting PCF value for radiusOfInterest
    pcf_value_at_radius = None
    if radiusOfInterest in r:
        index = np.where(r== radiusOfInterest)[0]
        if len(index) > 0:  
            pcf_value_at_radius = pcf[index][0]

    return pc, pcf_value_at_radius

def TCM(C_1, C_2, radiusOfInterest, pc, output_path):
    '''
    Function to evaluate and plot TCM 
    Args:
        C_1 : string
        C_2 : string
        radiusOfInterest : int
        pc : pointcloud
        output_path : string

    '''
    #computation of topographical correlation map 
    tcm = topographicalCorrelationMap(pc,'Celltype',C_1,'Celltype',C_2,radiusOfInterest,maxCorrelationThreshold=5.0,kernelRadius=150,kernelSigma=50,visualiseStages=False)
    #masked_tcm = np.ma.masked_where((tcm > -0.05) & (tcm < 0.05), tcm)

    #plot of the TCM
    plt.figure(figsize=(24,18), facecolor='none')
    l = int(np.ceil(np.max(np.abs([tcm.min(),tcm.max()]))))
    plt.imshow(tcm,cmap='RdBu_r',vmin=-l,vmax=l,origin='lower')
    plt.colorbar(label=f'$\Gamma_{{C_1}}$ $\Gamma_{{C_2}}$(r={radiusOfInterest})')
    ax = plt.gca()
    ax.grid(False)

    plt.tight_layout()
    plt.xticks([], [])
    plt.yticks([], [])

    plt.savefig(os.path.join(output_path, f"TCM.tif"), dpi=300, format="tiff", bbox_inches='tight')
    plt.close()

def append_to_csv(output_csv, codice_pz, roi_name, pcf_value_at_radius, count_C1, count_C2):
    '''
    Function to fill summary tsv
    Args:
        output_csv : string
        codice_pz : string
        numero_ROI : string
        pcf_value_at_radius : float
        count_C1 : string
        count_C2 : string
    '''
    with open(output_csv, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        pcf_value = pcf_value_at_radius[1].item() 
        row = [codice_pz, roi_name, pcf_value, count_C1, count_C2]
        writer.writerow(row)


### Algorithm for mathematical PCF evaluation
# Author: Bull, Joshua A., et al
# Work: "Extended correlation functions for spatial analysis of multiplex imaging data." Biological Imaging 4 (2024): e2.
# Git Repository: https://github.com/JABull1066/ExtendedCorrelationFunctions

class pointcloud:
    def __init__(self, 
    name, 
    points,
    domain=None, 
    unitOfLength=None
    ):
        self.name = name
        self.points = points #todo ensure points are an (n,d) numpy array for d = 2 or 3
        self.nPoints = np.shape(points)[0]
        self.dimension = np.shape(points)[1]       
        self.labels = {}
        self.nLabels = 0
        self.nLabels_categorical = 0
        self.nLabels_continuous = 0
        
        self.summaryStatistics = None

        if domain is None:
            # We estimate domain size by taking most extreme values and rounding up to nearest 1 unit
            maxDomain = np.ceil(np.max(self.points,axis=0))
            minDomain = np.floor(np.min(self.points,axis=0))
            self.domain = np.stack((minDomain,maxDomain)).transpose()
        else:
            self.domain = np.asarray(domain)
            if np.shape(self.domain)[0] != self.dimension:
                raise RuntimeError('Specified domain should be nDimensions by 2, specifying domain min and max values in each dimension')
            if np.shape(self.domain)[1] != 2:
                raise RuntimeError('Specified domain should be nDimensions by 2, specifying domain min and max values in each dimension')
            for d in range(self.dimension):
                if self.domain[d,0] >= self.domain[d,1]:
                    raise RuntimeError(f'Specified domain minimum value in dimension {d} ({self.domain[d,0]}) must be lower than maximum value ({self.domain[d,1]})')
        if self.dimension == 2:
            v = [[self.domain[0,0],self.domain[1,0]],
                 [self.domain[0,0],self.domain[1,1]],
                 [self.domain[0,1],self.domain[1,1]],
                 [self.domain[0,1],self.domain[1,0]]]
            self.boundaryPolygon = Polygon(v)
        self.domainVolume = np.prod(self.domain[:,1] - self.domain[:,0],axis=0)
        self.density = self.nPoints / self.domainVolume
        

    def __str__(self):
        return f"Name: {self.name}, nPoints: {self.nPoints}"

    def addLabels(self, labelName, labelType, labels,cmap=None):
        if self.nPoints != len(labels):
            raise ValueError(f"Expected a list of {self.nPoints} labels, received {len(labels)}")
        if labelType in ['categorical','continuous']:
            self.labels[labelName] = {'Type':labelType,
                                        'labels':labels}
            self.nLabels = self.nLabels + 1
            
            if labelType == 'categorical':
                self.nLabels_categorical = self.nLabels_categorical + 1
                unique = np.unique(labels)
                labelToInteger = {unique[v] : v for v in range(len(unique))}
                self.labels[labelName]['categories'] = unique
                self.labels[labelName]['nCategories'] = len(unique)
                self.labels[labelName]['labelToInteger'] = labelToInteger
                self.labels[labelName]['integerToLabel'] = {labelToInteger[v] : v for v in labelToInteger.keys()}
                self.labels[labelName]['numericalLabels'] = np.asarray([labelToInteger[labels[v]] for v in range(len(labels))])
                if cmap is None:
                    # Use default colormap
                    cmap = 'tab10'
                self.labels[labelName]['integerToColor'] = {v: plt.cm.tab20(v) for v in self.labels[labelName]['integerToLabel'].keys()}
                colArray = np.asarray([v for v in self.labels[labelName]['integerToColor'].values()])
                self.labels[labelName]['cmap'] = ListedColormap(colArray)
                
            else:
                self.nLabels_continuous = self.nLabels_continuous + 1
                self.labels[labelName]['numericalLabels'] = np.asarray(labels)
                self.labels[labelName]['cmap'] = 'plasma'
        else:
            raise ValueError('labelType must be categorical or continuous')

    def changeIndividualLabelColor(self, labelName, labelToUpdate, newColor):
        assert(len(newColor) == 4)
        labelIntegerValue = self.labels[labelName]['labelToInteger'][labelToUpdate]
        self.labels[labelName]['integerToColor'][labelIntegerValue] = newColor
        # Update cmap
        colArray = np.asarray([v for v in self.labels[labelName]['integerToColor'].values()])
        self.labels[labelName]['cmap'] = ListedColormap(colArray)

def generatePointCloud(name, points,domain=None,unitOfLength=None):
    pc = pointcloud(name, points, domain, unitOfLength)
    return pc

def visualisePointCloud(pc,labelForVisualisation=None,cmap=None,markerSize=None,vmin=None,vmax=None):
    from matplotlib import colors
    if pc.dimension != 2:
        raise RuntimeError('Visualisation currently only possible in 2D')
    if labelForVisualisation not in pc.labels.keys() and labelForVisualisation != None:
        raise ValueError('labelForVisualisation must be a label!')

    shuffleOrder = np.arange(len(pc.points))
    random.shuffle(shuffleOrder)

    if markerSize is None:
        markerSize = 20
    else:
        if isinstance(markerSize, (int,float)):
            if markerSize < 0:
                raise ValueError('markerSize must be a positive number')

    fig, ax = plt.subplots(figsize=(24,18))
    #plt.figure(figsize=(24,18))
    if labelForVisualisation is None:
        plt.scatter(pc.points[shuffleOrder,0],pc.points[shuffleOrder,1],s=markerSize,cmap=cmap)
    else:
        norm = None
        if cmap is None:
            cmap = pc.labels[labelForVisualisation]['cmap']
        labelType = pc.labels[labelForVisualisation]['Type']
        if labelType == 'categorical':
            nCategories = pc.labels[labelForVisualisation]['nCategories']
            # if cmap is None:
            #     cmap = pc.labels[labelForVisualisation]['cmap']
            cmap = plt.cm.get_cmap(cmap,nCategories)
            norm = colors.BoundaryNorm(np.arange(-0.5, nCategories+0.5, 1), cmap.N)
        plt.scatter(pc.points[shuffleOrder,0],pc.points[shuffleOrder,1],c=pc.labels[labelForVisualisation]['numericalLabels'][shuffleOrder],s=markerSize,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    plt.gca().axis('equal')
    plt.xlim(pc.domain[0])
    plt.ylim(pc.domain[1])
    if labelForVisualisation != None:
        cbar=plt.colorbar()#(label=labelForVisualisation)
        if labelType == 'categorical':
            cbar.set_ticks(list(pc.labels[labelForVisualisation]['labelToInteger'].values()))
            cbar.set_ticklabels(list(pc.labels[labelForVisualisation]['labelToInteger'].keys()))
    plt.tight_layout()
    plt.xticks([], [])
    plt.yticks([], [])

    cbar.ax.tick_params(labelsize=20)

    return plt.gcf(), plt.gca()

def returnIntersectionPoints(x0, y0, r, domainX, domainY):
    # Calculate the points of intersection between a circle of radius r centred at (x0,y0)
    # and the box boundaries x = domainX[0], y = domainY[0], x = domainX[1] and y = domainY[1]
    # This also includes corners which are within the domain

    # Find intersections with each of the 4 domain edges - gives max of 8 intersections. Then take floor/ceil to take points outside of domain to domain corners
    # We assume for now that domain edges are parallel to coordinate axes
    # i.e, vertices are [ [domainX[0],domainY[0]], [domainX[1],domainY[0]], [domainX[0],domainY[1]], [domainX[1],domainY[1]]  ]

    # Line of form ax + by = c
    # circle of form (x - x0)^2 + (y - y0)^2 = r^2
    # Need r^2 (a^2 + b^2) - (c - ax0 - by0)^2 > 0
    # See page 17 of https://www2.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
    def calculateUsefulValues(rSquared, a, b, c, x0, y0):
        cPrime = c - a*x0 - b*y0
        return rSquared*(a**2 + b**2) - cPrime**2, cPrime
    def getIntersections(val, a, b, cPrime, x0, y0):
        temp = np.sqrt(val)/(a**2 + b**2)
        point1 = [x0 + a*cPrime + b*temp, y0 + b*cPrime - a*temp]
        point2 = [x0 + a*cPrime - b*temp, y0 + b*cPrime + a*temp]
        return [point1, point2]

    rSquared = r**2
    intersectionPoints = []
    # Move around anti-clockwise from upper left corner
    # LEFT HAND SIDE
    # x = domainX[0], i.e. a=1, b=0, c=domainX[0]
    val, cPrime = calculateUsefulValues(rSquared, 1, 0, domainX[0], x0, y0)
    if val > 0:
        intersections = getIntersections(val, 1, 0, cPrime, x0, y0)
        intersectionPoints.extend([intersections[1],intersections[0]])
    # BOTTOM
    # y = domainY[0], i.e. a=0, b=1, c=domainY[0]
    val, cPrime = calculateUsefulValues(rSquared, 0, 1, domainY[0], x0, y0)
    if val > 0:
        intersections = getIntersections(val, 0, 1, cPrime, x0, y0)
        intersectionPoints.extend([intersections[1],intersections[0]])
    # RIGHT HAND SIDE
    # x = domainX[1], i.e. a=1, b=0, c=domainX[1]
    val, cPrime = calculateUsefulValues(rSquared, 1, 0, domainX[1], x0, y0)
    if val > 0:
        intersections = getIntersections(val, 1, 0, cPrime, x0, y0)
        intersectionPoints.extend(intersections)
    # BOTTOM
    # y = domainY[0], i.e. a=0, b=1, c=domainY[1]
    val, cPrime = calculateUsefulValues(rSquared, 0, 1, domainY[1], x0, y0)
    if val > 0:
        intersections = getIntersections(val, 0, 1, cPrime, x0, y0)
        intersectionPoints.extend(intersections)
    
    temp = [[np.max([np.min([point[0],domainX[1]]),domainX[0]]), np.max([np.min([point[1],domainY[1]]),domainY[0]])] for point in intersectionPoints]
    intersectionPoints = [temp[v] for v in range(len(temp)) if temp[v] not in temp[:v]]
    return intersectionPoints

def returnAreaOfCircleInDomain(x0, y0, r, domainX, domainY):
    intersectionPoints = returnIntersectionPoints(x0, y0, r, domainX, domainY)
    if not intersectionPoints:
        area = np.pi * r ** 2
    else:
        # Need to calculate area from intersection Points
        intersectionPoints.append(intersectionPoints[0])
        area = 0
        for v in range(len(intersectionPoints) - 1):
            a = intersectionPoints[v]
            b = intersectionPoints[v + 1]
            # Find out if this is a segment or a triangle
            isTriangle = False

            # Check if point b is anticlockwise from point a on the same line
            if a[0] == b[0] or a[1] == b[1]:
                if a[0] == b[0]:  # On a vertical line
                    if a[0] == domainX[0]:
                        # LHS
                        if b[1] < a[1]:
                            isTriangle = True
                    else:
                        # RHS
                        if b[1] > a[1]:
                            isTriangle = True
                else:  # On a horizontal line
                    if a[1] == domainY[0]:
                        # bottom
                        if b[0] > a[0]:
                            isTriangle = True
                    else:
                        # top
                        if a[0] > b[0]:
                            isTriangle = True

            # If points are on the same line moving anticlockwise, then return the area of the triangle formed by a, b and the centre
            if isTriangle:
                # Points share a border: return area of triangle between them
                area = area + 0.5 * np.abs(a[0] * (b[1] - y0) + b[0] * (y0 - a[1]) + x0 * (a[1] - b[1]))
            else:
                # Else, return the area of the circle segment between them
                # We need to be careful to take the angle between v1 and v2 in an anticlockwise direction
                v1 = [x0 - a[0], y0 - a[1]]
                v2 = [x0 - b[0], y0 - b[1]]

                theta = np.arctan2(v2[1], v2[0]) - np.arctan2(v1[1], v1[0])
                # Normalise to 0, 2pi
                if theta < 0:
                    theta = theta + 2 * np.pi

                area = area + 0.5 * theta * r ** 2
    return area

def returnAreaOfCircleInDomainAroundPoint(index, points, r, domainX, domainY):
    point = points[index,:]
    area = returnAreaOfCircleInDomain(point[0], point[1], r, domainX, domainY)
    return area

def plotTopographicalCorrelationMap(pc,topographicalCorrelationMap,ax=None,cmap='RdBu_r',colorbarLimit=None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if ax == None:
        plt.figure(figsize=(18,18), facecolor='none')
        ax = plt.gca()
    if colorbarLimit == None:
        colorbarLimit = int(np.ceil(np.max([topographicalCorrelationMap.min(),topographicalCorrelationMap.max()])))
    
    transparent_mask = np.ma.masked_where((topographicalCorrelationMap > -1) & (topographicalCorrelationMap<1), topographicalCorrelationMap)
    print(transparent_mask)


    extent = [pc.domain[0,0],pc.domain[0,1],pc.domain[1,0],pc.domain[1,1]]
    im = ax.imshow(topographicalCorrelationMap,origin='lower',cmap=cmap,extent=extent,vmin=-colorbarLimit,vmax=colorbarLimit, alpha=transparent_mask)
    
    im.patch.set_alpha(0)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.gcf().colorbar(im, cax=cax, orientation='vertical')
    return plt.gcf(), plt.gca()

def topographicalCorrelationMap(pc,labelNameA,labelA,labelNameB,labelB,radiusOfInterest,maxCorrelationThreshold=5.0,kernelRadius=150,kernelSigma=50,visualiseStages=False):
    
    for labelName in [labelNameA,labelNameB]:
        labelType = pc.labels[labelName]['Type']
        if labelType != 'categorical':
            raise RuntimeError(f'The label {labelName} is not a categorical label.')
    
    if labelA not in pc.labels[labelNameA]['categories']:
        raise RuntimeError(f'The category {labelA} is not associated with the label {labelNameA}.')
    if labelB not in pc.labels[labelNameB]['categories']:
        raise RuntimeError(f'The category {labelB} is not associated with the label {labelNameB}.')

    i_A = pc.labels[labelNameA]['labelToInteger'][labelA]
    i_B = pc.labels[labelNameB]['labelToInteger'][labelB]

    # Points to include A: All points within pc.domain
    # Points to include B: All points within pc.domain
    p_A = pc.points[pc.labels[labelNameA]['numericalLabels'] == i_A,:]
    p_B = pc.points[pc.labels[labelNameB]['numericalLabels'] == i_B,:]

    # Get areas around A, calculate pairwise A-B distances
    areas = []
    for i in range(len(p_A)):
        area = returnAreaOfCircleInDomainAroundPoint(i,p_A,radiusOfInterest,pc.domain[0],pc.domain[1])
        areas.append(area)
    density_B = np.shape(p_B)[0]/pc.domainVolume
    areas_A = np.asarray(areas)
    distances_AtoB = cdist(p_A, p_B, metric='euclidean')
    contributions = distances_AtoB <= radiusOfInterest
    
    BnearA_observed = np.sum(contributions,axis=1)/areas_A # observed per unit area
    marks = BnearA_observed/density_B
    
    
    if visualiseStages:
        s=100
        plt.figure(figsize=(20,20))
        plt.scatter(p_A[:,0],p_A[:,1],c=marks,cmap='viridis',s=s)
        plt.colorbar()
        plt.gca().axis('equal')
        
    # Map PCF interpretation to [-1,1]
    minCorrelationThreshold = 1/maxCorrelationThreshold
    
    transformedMarks = np.copy(marks)
    transformedMarks[transformedMarks < minCorrelationThreshold] = minCorrelationThreshold
    transformedMarks[transformedMarks > maxCorrelationThreshold] = maxCorrelationThreshold
    
    transformedMarks[transformedMarks<1] = -1/transformedMarks[transformedMarks<1]
    # That gives us values in [-maxPCFthreshold,-1] U [1,maxPCFthreshold]
    # Now map to [-1,1]
    transformedMarks[transformedMarks<0] = (transformedMarks[transformedMarks<0]+1)/(maxCorrelationThreshold-1)
    transformedMarks[transformedMarks>0] = (transformedMarks[transformedMarks>0]-1)/(maxCorrelationThreshold-1)
    
    if visualiseStages:
        plt.figure(figsize=(20,20))
        plt.scatter(p_A[:,0],p_A[:,1],c=transformedMarks,cmap='RdBu_r',vmin=-1,vmax=1,s=s)
        plt.colorbar()
        plt.gca().axis('equal')
                          
                        
    x, y = np.meshgrid(np.arange(-kernelRadius,kernelRadius+0.1,1),np.arange(-kernelRadius,kernelRadius+0.1,1))
    dst = np.sqrt(x*x + y*y)
    kernel = np.exp(-( dst**2 / ( 2.0 * kernelSigma**2 ) ) )
    
    xrange = [int(pc.domain[0][0])-kernelRadius, int(pc.domain[0][1])+1+kernelRadius]
    yrange = [int(pc.domain[1][0])-kernelRadius, int(pc.domain[1][1])+1+kernelRadius]
    heatmap = np.zeros(shape=(xrange[1]-xrange[0],yrange[1]-yrange[0]))
    
    def addWeightedContribution(heatmap, weight, coordinate, xrange, yrange, kernel,kernelRadius):
        x0 = int(coordinate[0]) - kernelRadius - xrange[0]
        x1 = x0 + 2*kernelRadius + 1
        y0 = int(coordinate[1]) - kernelRadius - yrange[0]
        y1 = y0 + 2*kernelRadius + 1
        heatmap[x0:x1,y0:y1] = heatmap[x0:x1,y0:y1] + kernel*weight
        return heatmap
    
    for i in range(len(p_A)):
        coordinate = p_A[i,:]
        weight = transformedMarks[i]
        heatmap = addWeightedContribution(heatmap, weight, coordinate, xrange, yrange, kernel, kernelRadius)
    
    topographicalCorrelationMap = heatmap[kernelRadius:-kernelRadius,kernelRadius:-kernelRadius]

    if visualiseStages:
        l = int(np.ceil(np.max([topographicalCorrelationMap.min(),topographicalCorrelationMap.max()])))
        fig, ax = plotTopographicalCorrelationMap(pc,topographicalCorrelationMap.T,ax=None,cmap='RdBu_r',colorbarLimit=l)
    
    return topographicalCorrelationMap.T

def crossPCF(distances_AtoB, areas_A, density_B, maxR, annulusStep, annulusWidth):
    N_A = np.shape(distances_AtoB)[0]

    PCF_radii_lower = np.arange(0, maxR + annulusStep, annulusStep)
    PCF_radii_upper = np.arange(annulusWidth, maxR + annulusStep + annulusWidth, annulusStep)

    crossPCF_AtoB = np.ones(shape=(len(PCF_radii_lower),1))
    contributions = np.zeros(shape=(N_A,len(PCF_radii_lower)))
    for annulus in range(len(PCF_radii_lower)):
        inner = PCF_radii_lower[annulus]
        outer = PCF_radii_upper[annulus]

        # Find pairwise distances within this radius
        distanceMask = np.logical_and((distances_AtoB > inner),(distances_AtoB <= outer))
        for i in range(N_A):
            # For each point in pA
            # Find pairwise distances to points in pB within this radius
            fillIndices = np.where(distanceMask[i,:])[0]
            contribution = len(fillIndices)/(density_B*areas_A[i,annulus])
            crossPCF_AtoB[annulus] = crossPCF_AtoB[annulus] + contribution
            contributions[i,annulus] = contributions[i,annulus] + contribution
        crossPCF_AtoB[annulus] = crossPCF_AtoB[annulus] / N_A
    return PCF_radii_lower, crossPCF_AtoB, contributions

def pairCorrelationFunction(pc,labelName,categoriesToPlot,maxR=0.5,annulusStep=0.025,annulusWidth=0.025):
# First we check that the chosen label is categorical
    labelType = pc.labels[labelName]['Type']
    if labelType != 'categorical':
        raise RuntimeError(f'The label {labelName} is not a categorical label.')
    categories = pc.labels[labelName]['categories']
    labelA = categoriesToPlot[0]
    labelB = categoriesToPlot[1]
    if labelA not in categories:
        raise RuntimeError(f'The category {labelA} is not associated with the label {labelName}.')
    if labelB not in categories:
        raise RuntimeError(f'The category {labelB} is not associated with the label {labelName}.')

    
    i_A = pc.labels[labelName]['labelToInteger'][labelA]
    i_B = pc.labels[labelName]['labelToInteger'][labelB]
    
    # Points to include A: All points within pc.domain
    # Points to include B: All points within pc.domain
    p_A = pc.points[pc.labels[labelName]['numericalLabels'] == i_A,:]
    p_B = pc.points[pc.labels[labelName]['numericalLabels'] == i_B,:]
    if np.shape(p_A)[0] == 0:
        raise RuntimeError(f'No cells with {labelA} found within PCF domain')
    if np.shape(p_B)[0] == 0:
        raise RuntimeError(f'No cells with {labelB} found within PCF domain')
    # Get annulus areas (within domain) around p_A
    areas_A = getAnnulusAreasAroundPoints(p_A, maxR, annulusStep,annulusWidth, pc.domain)
    density_B = np.shape(p_B)[0]/pc.domainVolume
    
    distances_AtoB = cdist(p_A, p_B, metric='euclidean')
    radii, g, contributions = crossPCF(distances_AtoB, areas_A, density_B, maxR, annulusStep, annulusWidth)

    return radii, g, contributions

def getAnnulusAreasAroundPoints(points_i, maxR, annulusStep, annulusWidth, domain):
    # We want to populate a table the same size as distances, which contains the area of the annulus containing that contribution
    # i.e., "at distance D(i->j) from point i, what is area of containing annulus?"
    domainX = domain[0,:]
    domainY = domain[1,:]
    vfunc_returnAreaOfCircleInDomainAroundPoint = np.vectorize(returnAreaOfCircleInDomainAroundPoint,excluded=['points','domainX','domainY'])
    PCF_radii_lower = np.arange(0, maxR+annulusStep, annulusStep)
    PCF_radii_upper = np.arange(annulusWidth, maxR + annulusWidth + annulusStep, annulusStep)
    # PCF_radii_lower = np.arange(0, maxR, dr)
    # PCF_radii_upper = np.arange(dr, maxR + dr, dr)

    allAreas = np.zeros(shape=(len(points_i),len(PCF_radii_lower)))

    for annulus in range(len(PCF_radii_lower)):
        inner = PCF_radii_lower[annulus]
        outer = PCF_radii_upper[annulus]

        areas_in = vfunc_returnAreaOfCircleInDomainAroundPoint(index=np.arange(len(points_i)), points=points_i, r=inner, domainX=domainX, domainY=domainY)
        areas_out = vfunc_returnAreaOfCircleInDomainAroundPoint(index=np.arange(len(points_i)), points=points_i, r=outer, domainX=domainX, domainY=domainY)

        areas = areas_out - areas_in
        if not np.all(areas >= 0):
            raise RuntimeError(f'Negative areas calculated for point {np.argwhere(areas < 0)}.')
        allAreas[:,annulus] = areas
    return allAreas