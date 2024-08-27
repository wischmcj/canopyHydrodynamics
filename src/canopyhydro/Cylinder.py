"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from canopyhydro.configuration import qsm_cols
from canopyhydro.DataClasses import Projection
from canopyhydro.geometry import draw_cylinders_3D, draw_cyls, get_projection

log = logging.getLogger("model")

NAME = "Cylinder"


def create_cyl(arr: np.array):
    """Creates a cylinder based off of an input array.
        Enables user defined file structure via configuration yml file

    Args:
        arr (np.array):
            Array containing attribute data for a single cylinder.
            Attributes must be ordered as per qsm_cols in canopyhydro_config.yml
    Examples:
        >>> import numpy as np
        >>> from canopyhydro.configuration import qsm_cols
        >>> print(qsm_cols)
        {'cyl_id': 1, 'parent_id': 2, 'x': [3, 6], 'y': [4, 7], 'z': [5, 8],
        'radius': 9, 'volume': 10, 'length': 12, 'segment_id': 15, 'branch_order': 20,
        'reverse_branch_order': 21, 'branch_id': 24}
        >>> myCyl = create_cyl(
        >>>    np.array([np.nan, 1, 0, -2.372914, 2.875943, -0.55, -2.382034, 2.887813,
        >>>    -0.452896, 0.277545, 0.023777, 4.755711, 0.098251, 3880.839667,
        >>>    np.nan, 0, -1, 0.240459, 4.604282, 3880.04742, 0, 41, 84.320816,
        >>>    7110, 0, 0, np.nan, 0, 0, 0, 0.01, 0.401606, 0, 0.01, 0.401606])
        >>> )
        >>> print(myCyl)
        Cylinder( cyl_id=1.0, x=[-2.372914 -2.382034], y=[2.875943 2.887813], z=[-0.55     -0.452896], radius=0.277545, length=0.098251, branch_order=0.0, branch_id=0.0, volume=0.023777, parent_id=0.0, reverse_branch_order=41.0, segment_id=0.0
    """
    cols = qsm_cols
    print(cols)
    attrs = {k: arr[v] for (k, v) in cols.items()}
    cyl = Cylinder(**attrs)
    # cyl.initialize(arr, cols)
    return cyl


@dataclass
class Cylinder:
    """
    The Cylinder class is used to represent the 3-D cylinders that make up a QSM.
    Contains several wrappers for functions in 'geometry'.

    Attributes:
        cyl_id (int):
            The ID of the cylinder.
        x (np.ndarray[np.float32]):
            The x-coordinates of the cylinder vertices. Length 2, ordered (start, end)
        y (np.ndarray[np.float32]):
            The y-coordinates of the cylinder vertices. Length 2, ordered (start, end)
        z (np.ndarray[np.float32]):
            The z-coordinates of the cylinder vertices. Length 2, ordered (start, end)
        radius (np.float32):
            The radius of the cylinder.
        length (np.float32):
            The length of the cylinder.
        branch_order (int):
            The branch order of the cylinder. See QSM documentation for more information.
        branch_id (int):
            The ID of the branch the cylinder belongs to. See QSM documentation for more information.
        volume (np.float32):
            The volume of the cylinder.
        parent_id (int):
            The ID of the parent cylinder. -1 if cylinder is a root cylinder
        reverse_branch_order (int):
            The reverse branch order of the cylinder. See QSM documentation for more information.
        segment_id (int):
            The ID of the segment the cylinder belongs to. See QSM documentation for more information.
        projected_data (dict(Projection)):
            A dictionary containing projected data of the cylinder
            on the 'XY', 'XZ' and 'YZ' planes. Defaults to an empty dictionary.
        flow_id (int):
            The ID of the flow to which the cylinder contributes intercepted precipitation.
            Defaults to None.
        flow_type (str):
            The type of flow associated with the cylinder (drip or stem). Defaults to None.
        drip_node (int):
            The ID of the drip node associated with the cylinder.
            Defined by an ID correlating to the node_id of the node in the
            containing CylinderCollection's graph . Defaults to None.
        begins_at_drip_point (bool):
            Indicates if (x[0],y[0],z[0]) at a drip point. Defaults to None.
        begins_at_divide_point (bool):
            Indicates if (x[0],y[0],z[0]) at a dividepoint. Defaults to None.
        dx (np.float32):
            The change in x-coordinate: x[1] - x[0]. Defaults to 0.
        dy (np.float32):
            The change in y-coordinate: y[1] - y[0]. Defaults to 0.
        dz (np.float32):
            The change in z-coordinate: z[1] - z[0]. Defaults to 0.
        surface_area (np.float32):
            The surface area of the cylinder. Defaults to 0.0.
        sa_to_vol (np.float32):
            The ratio of surface area to volume of the cylinder. Defaults to 0.0.
        slope (np.float32):
            The slope of the cylinder, the rise over run in 3 dimensions. Defaults to 0.0.
        is_stem (bool):
            Indicates if the cylinder contributes to stem flow. Defaults to False.
    """

    cyl_id: int
    x: np.ndarray[np.float32]
    y: np.ndarray[np.float32]
    z: np.ndarray[np.float32]
    radius: np.float32
    length: np.float32
    branch_order: int
    branch_id: int
    volume: np.float32
    parent_id: int
    reverse_branch_order: int
    segment_id: int

    projected_data: dict[str, Projection] = field(default_factory=dict)
    flow_id: int = None
    flow_type: str = None
    drip_node: int = None
    begins_at_drip_point: bool = None
    begins_at_divide_point: bool = None

    dx: np.float32 = 0
    dy: np.float32 = 0
    dz: np.float32 = 0

    surface_area: np.float32 = 0.0
    sa_to_vol: np.float32 = 0.0
    slope: np.float32 = 0.0

    is_stem: bool = False

    def __repr__(self):
        """Defines a readable string to represent a cylinder object
        *Utilized in tests to compare expected and actual results
        """
        return f"Cylinder( cyl_id={self.cyl_id}, x={self.x}, y={self.y}, z={self.z}, radius={self.radius}, length={self.length}, branch_order={self.branch_order}, branch_id={self.branch_id}, volume={self.volume}, parent_id={self.parent_id}, reverse_branch_order={self.reverse_branch_order}, segment_id={self.segment_id}"

    def __eq__(self, other):
        """Defines the minimum requirements for equality between two cylinders.
        Allows use of '=' for comparing cylinders
        """
        return type(self) == type(other) and self.__repr__() == other.__repr__()

    def __post_init__(self):
        """Initializes the Cylinder object after all attributes are set"""
        self.initialize()

    def initialize(self):
        """Initializes remaining attributes based off of the
        attributes provided at object creation

        References:
            - self.x, self.y, self.z
            - self.dx, self.dy, self.dz
            - self.radius, self.length
            - self.volume
            - self.surface_area
            - self.sa_to_vol
            - self.angle
            - self.xy_area

        """

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        radius = self.radius
        length = self.length
        self.surface_area = (
            2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
        )

        self.sa_to_vol = 0 if self.volume == 0 else self.surface_area / self.volume
        run = np.sqrt(self.dx**2 + self.dy**2)
        self.angle = (
            np.arctan(0)
            if run == 0
            else np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
        )
        self.xy_area = 0

    def get_projection(self, plane="XY"):
        """Calculates the projection of the cylinder on the requested plane

        Args:
            self (Cylinder):
                Cylinder to be projected
            plane (str, optional):
                Plane on which to project the cylinder. Defaults to "XY".

        References:
            projected_data.polygon (shapely.Polygon):
              represents the shape formed by projecting
              the cylinder (self) onto the given plane
            projected_data.base_vector (list[float]):
              cylinder vector *after projection
              (oriented away from cylinder base)
            projected_data.anti_vector (list[float]):
              cylinder vector *after projection
              (oriented towards cylinder base)
            projected_data.angle (float):
              angle (in radians) of cylinder vector with
              given axis. Between -pi/2 and pi/2.
            projected_data.area (float) - area of polygon (see above)
            projected_data.xy_area: float
              - polygon formed by projecting the cylinder
                  onto the XY plane

        Returns:
            (shapely.Polygon): projected_data[polygon]

        Examples:
            >>> import numpy as np
            >>> from canopyhydro.Cylinder import Cylinder
            >>> cyl = Cylinder(1, np.array([0, 1]), np.array([0, 1]), np.array([0, 1]), 1, 1, 0, 0, 1, 0, 0, 0)
            >>> cyl.get_projection("XY")
            >>> print(cyl.projected_data['XY']['polygon'])
            POLYGON ((0.7055547485197222 -0.7086486608663443, 0.6927794439882502 -0.7206122246452861, 0.6793868045462438 -0.7318989317552199, 0.6654090505126061 -0.7424887326410364, 0.6508814456329607 -0.7523657839197849, 0.6358418765269868 -0.761518538618793, 0.6203303971002122 -0.7699397615971714, 0.6043887500439113 -0.777626472292711, 0.5880598778050475 -0.7845798192275877, 0.5713874351273125 -0.7908048926864933, 0.5544153144933113 -0.7963104835647774, 0.5371871946141021 -0.8011087975191887, 0.5197461206122052 -0.8052151342206525, 0.5021341228348492 -0.8086475417192888, 0.4843918794234269 -0.8114264557256392, 0.4665584259537658 -0.8135743330482637, 0.4486709137372394 -0.8151152875793871, 0.4307644168049181 -0.8160747361660723, 0.4128717862361091 -0.81647906052316, 0.3950235493693495 -0.8163552901089949, 0.3772478505606316 -0.8157308096593167, 0.3595704295270235 -0.8146330939104094, 0.3420146329184335 -0.813089470978599, 0.3246014545714764 -0.811126914925819, 0.3073495998873245 -0.8087718672450097, 0.290275569907372 -0.8060500863494255, 0.2733937609038954 -0.8029865236433437, 0.2567165756266195 -0.7996052243790617, 0.2402545427223602 -0.7959292512531837, 0.2240164412495102 -0.7919806285483623, 0.2080094276222298 -0.7877803045682026, 0.1922391627249603 -0.7833481301263573, 0.1767099373244321 -0.7787028509202096, 0.1614247942654083 -0.7738621117306963, 0.1463856462625561 -0.7688424705301826, 0.1315933883913461 -0.7636594207391983, 0.1170480046348049 -0.7583274200414502, 0.1027486680610392 -0.7528599243378058, 0.0886938343903684 -0.7472694255884285, 0.0748813288632188 -0.7415674924539188, 0.061308426443543 -0.7357648127983645, 0.0479719254905637 -0.7298712372577524, 0.0348682151073508 -0.7238958232052886, 0.0219933364311366 -0.7178468785603557, 0.0093430381703601 -0.7117320049902478, -0.0030871732801796 -0.7055581401438671, -0.0153019887994878 -0.6993315986349682, -0.0273062562232358 -0.6930581115600766, -0.0391049443978219 -0.6867428643938861, -0.0507031112609087 -0.680390533153648, -0.0621058756146832 -0.674005318764837, -0.0733183922721946 -0.6675909795940658, -0.0843458302733461 -0.661150862142771, -0.0951933538847147 -0.6546879299173791, -0.1058661061157358 -0.6482047905092289, -0.1163691945023672 -0.641703720931191, -0.1267076789278345 -0.6351866912682232, -0.136886561268094 -0.6286553867066296, -0.1469107766670601 -0.6221112280119577, -0.1567851862632475 -0.6155553905287505, -0.1665145712052036 -0.6089888217770355, -0.1761036278078615 -0.6024122577208589, -0.1855569637157113 -0.5958262377835759, -0.1948790949514776 -0.5892311186832065, -0.2040744437407678 -0.5826270871591742, -0.213147337014027 -0.5760141716592699, -0.2221020054970459 -0.5693922530529014, -0.230942583310351 -0.5627610744336571, -0.2396731080060256 -0.5561202500710708, -0.2482975209780038 -0.5494692735682325, -0.2568196681886216 -0.5428075252786534, -0.265243301160299 -0.5361342790325831, -0.2735720781867121 -0.5294487082198208, -0.2818095657226979 -0.5227498912730146, -0.2899592399165438 -0.5160368165924882, -0.2980244882522005 -0.5093083869508169, -0.3060086112724441 -0.5025634234126783, -0.3139148243570635 -0.4958006688029721, -0.3217462595328715 -0.4890187907537801, -0.3295059672946971 -0.4822163843585018, -0.337196918418591 -0.475391974459369, -0.344822005750257 -0.468544017592593, -0.3523840459532637 -0.4616709036135613, -0.359885781202887 -0.4547709570228152, -0.3673298808125219 -0.4478424380119915, -0.3747189427804962 -0.4408835432474799, -0.3820554952458164 -0.4338924064082732, -0.3893419978419239 -0.4268670984933031, -0.3965808429379147 -0.4198056279125201, -0.4037743567569047 -0.4127059403750534, -0.4109248003613197 -0.4055659185869789, -0.4180343704948377 -0.3983833817705406, -0.4251052002705446 -0.3911560850161054, -0.4321393596945685 -0.383881718477677, -0.4391388560140332 -0.3765579064224782, -0.4461056338776438 -0.3691822061448938, -0.4530415752965515 -0.3617521067550056, -0.4599484993923817 -0.3542650278519916, -0.4668281619184035 -0.3467183180928712, -0.4736822545388206 -0.339109253667407, -0.4805124038500116 -0.3314350366904837, -0.4873201701263032 -0.3236927935239349, -0.4941070457714587 -0.3158795730406438, -0.5008744534555515 -0.3079923448447743, -0.5076237439152422 -0.3000279974632291, -0.5143561933936759 -0.2919833365249203, -0.5210730006942992 -0.2838550829461442, -0.5277752838208069 -0.2756398711423602, -0.53446407617322 -0.2673342472889547, -0.5411403222677194 -0.2589346676562012, -0.5478048729453583 -0.2504374970465884, -0.5544584800321164 -0.2418390073660721, -0.561101790409973 -0.2331353763645797, -0.5677353394557562 -0.2243226865853769, -0.5743595438014955 -0.2153969245676656, -0.580974693366879 -0.2063539803521097, -0.5875809426112121 -0.1971896473449251, -0.5941783009490393 -0.1878996226027388, -0.6007666222703562 -0.1784795076077271, -0.6073455935031509 -0.1689248096105697, -0.6139147221529605 -0.1592309436276274, -0.6204733227512649 -0.149393235188444, -0.6270205021419827 -0.1394069239403035, -0.6335551435332125 -0.1292671682281142, -0.6400758892397821 -0.1189690507804404, -0.6465811220413709 -0.1085075856460081, -0.6530689450811008 -0.0978777265395345, -0.6595371602308326 -0.0870743767712113, -0.6659832448522456 -0.0760924009505614, -0.6724043268874385 -0.0649266386726187, -0.6787971582196869 -0.0535719204122682, -0.6851580862545552 -0.0420230858709856, -0.691483023684323 -0.030275005038799, -0.6977674164152208 -0.018322602252764, -0.7040062096579403 -0.0061608835510669, -0.7101938122080033 0.0062150323614731, -0.7163240589746696 0.0188098792060368, -0.7223901718559587 0.0316282039382399, -0.7283847191040086 0.0446743229220967, -0.7342995733803309 0.057952273116563, -0.7401258687655152 0.0714657579257982, -0.745853957063534 0.0852180873770077, -0.7514733638278601 0.0992121123130888, -0.7569727446358806 0.1134501523237356, -0.7623398422500517 0.1279339171907024, -0.7675614454291246 0.1426644216935309, -0.7726233502902727 0.1576418937141919, -0.777510325272263 0.1728656756960196, -0.7822060809093063 0.1883341196571653, -0.7866932457923251 0.2040444761345575, -0.7909533502654222 0.2199927776437158, -0.794966819575187 0.2361737174846068, -0.7987129783525699 0.252580525005119, -0.8021700684529729 0.2692048387511439, -0.80531528229976 0.2860365792834941, -0.8081248139574913 0.3030638238222093, -0.8105739301899368 0.3202726852807528, -0.8126370637189466 0.3376471986651436, -0.814287930777089 0.3551692182215297, -0.8154996748228801 0.3728183291011376, -0.8162450379463569 0.3905717776508622, -0.8164965610207151 0.4084044247040278, -0.8162268130423783 0.4262887264086602, -0.8154086493423686 0.4441947471576499, -0.8140154974488681 0.4620902090435216, -0.8120216683465209 0.4799405819193441, -0.8094026897356591 0.4977092175803609, -0.8061356566799797 0.515357530770148, -0.8021995937928458 0.532845228653828, -0.7975758219111152 0.5501305890972729, -0.7922483211122132 0.5671707865712363, -0.786204081022744 0.5839222628081372, -0.7794334287247852 0.6003411375420037, -0.7719303242641451 0.616383652842104, -0.763692613866814 0.632006642806133, -0.7547222315203778 0.6471680188162006, -0.7450253405956502 0.6618272592883419, -0.7346124086591418 0.6759458919638907, -0.7234982105148698 0.689487956381826, -0.7117017567372352 0.7024204342918837, -0.7070306424686519 0.7070306424686519, -0.7071067811865476 0.7071067811865476, 0.2928932188134524 1.7071067811865475, 0.292969357531348 1.7070306424686519, 0.3007538525912017 1.714713636442557, 0.313841646670986 1.7263415353923266, 0.3275310705289108 1.7372820356944092, 0.3417883322084627 1.747517174907592, 0.356576801490866 1.7570332512720226, 0.371857449443115 1.7658208764237213, 0.3875893170067164 1.773874954061725, 0.4037300003051736 1.7811945878880078, 0.4202361403636406 1.7877829242886305, 0.4370639054681907 1.7936469370153314, 0.4541694553806597 1.7987971624907408, 0.471509377979925 1.8032473952616899, 0.489041090518847 1.8070143535605634, 0.5067231994593613 1.8101173249305886, 0.5245158146712108 1.8125778014774983, 0.5423808155576559 1.8144191135941967, 0.5602820683252652 1.8156660700437288, 0.578185595084233 1.816344611158542, 0.5960596967104048 1.8164814806971856, 0.613875032398947 1.8161039206619238, 0.6316046595883158 1.8152393921800862, 0.6492240384415988 1.8139153244337884, 0.6667110053609246 1.8121588926194068, 0.6840457201070964 1.8099968250503737, 0.7012105910325372 1.807455238793954, 0.71819018274397 1.804559502655354, 0.7349711102240057 1.801334125884361, 0.75154192308756 1.7978026706693055, 0.7678929832558364 1.7939876862858326, 0.7840163389195337 1.7899106626674217, 0.7999055972517952 1.7855920011440585, 0.8155557979340072 1.7810510001387028, 0.8309632891839002 1.776305853703123, 0.8461256076323035 1.7713736609018986, 0.8610413630863489 1.7662704442042227, 0.875710128944702 1.761011175207885, 0.8901323387945563 1.7556098061905354, 0.9043091895193891 1.7500793061537925, 0.9182425510787038 1.7444317001912686, 0.9319348829834075 1.738678111168737, 0.9453891573800489 1.732828802851122, 0.9586087885706603 1.7268932237454224, 0.9715975687292824 1.7208800510503117, 0.984359609528348 1.7147972342119306, 0.9968992893551921 1.7086520376815417, 1.009221205778446 1.7024510825548362, 1.0213301329137072 1.6962003868455502, 1.03323098333563 1.689905404208528, 1.044928774187713 1.6835710609804537, 1.056428597150096 1.6772017914510826, 1.0677355919383413 1.670801571314928, 1.0788549230214501 1.6643739492838692, 1.0897917592643616 1.6579220768659142, 1.1005512562182267 1.651448736335148, 1.1111385408002799 1.6449563669334186, 1.121558698123705 1.6384470893562417, 1.1318167602561795 1.631922728584236, 1.1419176967035314 1.6253848351277167, 1.151866406431961 1.61883470475623, 1.1616677112584637 1.6122733967872591, 1.1713263504543203 1.6057017510093408, 1.1808469764208023 1.599120403314707, 1.1902341513094987 1.592529800115559, 1.1994923444719707 1.5859302116163443, 1.208625930634746 1.5793217440121707, 1.2176391887060567 1.57270435068084, 1.2265363011302155 1.5660778424330701, 1.2353213537141654 1.5594418968823762, 1.2439983358585964 1.5527960669928793, 1.2525711411331217 1.5461397888600703, 1.2610435681414325 1.5394723887763258, 1.2694193216281078 1.5327930896297834, 1.2777020137839628 1.5261010166820874, 1.2858951657114353 1.5193952027674997, 1.2940022090156669 1.5126745929529952, 1.3020264874906133 1.5059380486961969, 1.3099712588727794 1.4991843515353898, 1.3178396966380659 1.492412206343385, 1.325634891819738 1.4856202441746607, 1.3333598548277552 1.4788070247330372, 1.3410175172516086 1.4719710384850944, 1.348610733630488 1.4651107084426507, 1.3561422831760046 1.4582243916358582, 1.3636148714338916 1.4513103802968557, 1.3710311318720891 1.444366902772434, 1.378393627383419 1.4373921241828063, 1.3857048516916757 1.4303841468423566, 1.3929672306504137 1.4233410104571256, 1.400183123424021 1.4162606921128182, 1.4073548235408244 1.4091411060662453, 1.4144845598079971 1.4019801033523769, 1.4215744970779283 1.394775471218551, 1.4286267368554835 1.3875249323968777, 1.4356433177352237 1.3802261442254897, 1.4426262156571752 1.3728766976290254, 1.4495773439691424 1.36547411596859, 1.4564985532828494 1.358015853771428, 1.4633916311103494 1.3504992953506685, 1.4702583012662032 1.34292175332577, 1.477100223019839 1.335280467054708, 1.4839189899813228 1.3275726009895321, 1.4907161287024315 1.3197952429676691, 1.4974930969734763 1.3119454024522819, 1.5042512817947338 1.3040200087361404, 1.5109919969996237 1.2960159091248127, 1.517716480504905 1.287929867116582, 1.524425891161167 1.2797585605983501, 1.5311213051747354 1.2714985800789316, 1.5378037120698291 1.263146426983586, 1.5444740101573553 1.2546985120364398, 1.551133001474157 1.2461511537606063, 1.5577813861537941 1.2375005771293914, 1.5644197561871 1.228742912405994, 1.5710485885277579 1.219874194213621, 1.5776682374950797 1.2108903608829822, 1.584278926422992 1.2017872541297479, 1.5908807385010073 1.1925606191208094, 1.5974736067497277 1.1832061049951044, 1.604057303070199 1.1737192659124374, 1.610631426303307 1.1640955627121583, 1.6171953892324442 1.154330365272848, 1.623748404459941 1.1444189556743072, 1.630289469085405 1.1343565322742313, 1.6368173481122401 1.1241382148239967, 1.6433305565074043 1.1137590507610018, 1.6498273398391015 1.1032140228290355, 1.6563056534178071 1.0924980581931358, 1.6627631398680773 1.081606039231364, 1.6691971050622976 1.0705328162027226, 1.6756044923532536 1.0592732220080454, 1.6819818550505752 1.047822089278847, 1.6883253270971965 1.036174270047658, 1.694630591916554 1.0243246582719474, 1.7008928494198972 1.012268215501945, 1.7071067811865475 1, 1.7071067811865475 0.9999999999999998, 1.7132665138589451 0.9875151996348865, 1.7193655808296884 0.9748091688878777, 1.7253968823404577 0.9618774703174509, 1.7313526441635505 0.9487159208349948, 1.7372243750968257 0.9353206431433718, 1.7430028235730284 0.9216881226820279, 1.7486779337656866 0.9078152704045962, 1.7542388016668713 0.8936994916853315, 1.7596736317176722 0.8793387616068311, 1.7649696946906661 0.8647317068205472, 1.770113287654922 0.8498776940905957, 1.7750896969976253 0.8347769255272534, 1.7798831656310838 0.8194305403860189, 1.7844768656775423 0.8038407231481788, 1.7888528780938655 0.7880108174065364, 1.7929921808693186 0.7719454448530824, 1.7968746475965824 0.7556506284022838, 1.8004790583713068 0.739133918184136, 1.8037831251096352 0.722404518806347, 1.806763533475093 0.7054734159185241, 1.8093960036629888 0.6883534997189222, 1.811655372287373 0.6710596826354873, 1.8135156975366242 0.653609008000322, 1.81495038959218 0.6360207461371238, 1.8159323680242552 0.6183164739156477, 1.8164342474734188 0.6005201335208182, 1.8164285523851826 0.5826580659656443, 1.8158879608785283 0.564759014778021, 1.8147855769971704 0.5468540953447967, 1.8130952296215144 0.5289767256342404, 1.810791795227253 0.5111625144692783, 1.8078515404927495 0.4934491042111024, 1.804252479523993 0.4758759656487898, 1.799974739237794 0.4584841440751037, 1.795000925287278 0.4413159569448355, 1.7893164799034582 0.4244146451245806, 1.7829100222420806 0.4078239814963779, 1.7757736613440929 0.3915878424977104, 1.7679032717111878 0.3757497499756175, 1.7592987218199836 0.3603523923999496, 1.7499640466822795 0.3454371359115821, 1.7399075568092086 0.3310435367705419, 1.729141877626351 0.3172088674237017, 1.71768391545405 0.3039676685611681, 1.7070712498977356 0.2929287501022643, 1.7071067811865475 0.2928932188134524, 0.7071067811865476 -0.7071067811865476, 0.7070712498977356 -0.7070712498977356, 0.7055547485197222 -0.7086486608663443))
            >>> print(cyl.projected_data['XY']['base_vector'])
            [0.57735027 0.57735027 0.57735027]
            >>> print(cyl.projected_data['XY']['anti_vector'])
            [-0.57735027 -0.57735027 -0.57735027]
            >>> print(cyl.projected_data['XY']['angle'])
            0.6154797086703873
            >>> print(cyl.projected_data['XY']['area'])
            4.642087600836025
        """
        if plane == "XY":
            magnitude = [self.dx, self.dy, self.dz]
            ranges = [self.x, self.y, self.z]
        elif plane == "XZ":
            magnitude = [self.dx, self.dz, self.dy]
            ranges = [self.x, self.z, self.y]
        else:
            magnitude = [self.dy, self.dz, self.dx]
            ranges = [self.y, self.z, self.x]

        vector = list(map(np.transpose, ranges))

        projection = get_projection(vector, magnitude, self.radius)
        self.projected_data[plane] = projection
        if plane == "XY":
            self.xy_area = self.projected_data["XY"]["area"]
        return projection["polygon"]

    def draw(self, plane: str = "XY", **kwargs):
        """A wrapper around the draw_cyl function allowing
        more readable code for drawing a single cylinder.
        e.g. For some Cylinder 'cyl' cyl.draw() replaces geometry.draw_cyls([cyl])

        Args:
            plane:
                The projection of the cylinder to draw.
                'XY, 'XZ', or 'YZ'. Defaults to "XY".
        """
        poly = self.projected_data[plane]["polygon"]
        draw_cyls([poly], **kwargs)

    def draw_3D(self, **kwargs):
        """Draws the cylinder in 3D space"""
        vector_start_end = np.array(
            [
                np.array([self.x[0], self.y[0], self.z[0]]),
                np.array([self.x[1], self.y[1], self.z[1]]),
            ]
        )
        fig = draw_cylinders_3D([self.radius], [vector_start_end], **kwargs)
        return fig
