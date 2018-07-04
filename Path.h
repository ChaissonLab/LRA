#ifndef PATH_H_
#define PATH_H_

enum Arrow { Diagonal, Up, Left, 
						 AffineInsUp, AffineInsOpen, AffineInsClose,
						 AffineDelLeft, AffineDelOpen, AffineDelClose,
						 AffineHPInsUp, AffineHPInsOpen, AffineHPInsClose,
						 NoArrow,
						 DiagonalXYZ,
						 InsertX,InsertY,InsertZ, // imply diagonal yz/xz/xy
						 DiagonalXY, DiagonalYZ, DiagonalXZ,  // imply insertion of Z/X/Y,
             //
             // These are used to denote an affine gap has been closed
             // from a different matrix.  This is used in OneGap alignment.
             // 
             AffineLongDelLeft, AffineLongDelClose, 
             AffineLongIns, AffineLongInsClose,
             Star
};

enum MatrixLabel {Match, AffineHPIns, AffineIns, AffineDel, AffineHPDel};

#endif
