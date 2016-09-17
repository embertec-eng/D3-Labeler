(function() {
    d3.labeler = function() {
        var lab = [],
            anc = [],
            boundary = { // hard boundary
                left: 0,
                top: 0,
                right: 1,
                bottom: 1
            };
            labeler = {};

        var max_move = 5.0,
            max_angle = 0.5,
            acc = 0;
            rej = 0;

        // weights
        var w_len = 0.2, // leader line length
            w_inter = 1.0, // leader line intersection
            w_lab2 = 30.0, // label-label overlap
            w_lab_anc = 30.0; // label-anchor overlap
            w_orient = 3.0; // orientation bias

        // booleans for user defined functions
        var user_energy = false,
            user_schedule = false;

        var user_defined_energy,
            user_defined_schedule;

        function label_boundary(lbl) {
            return {
                left: lbl.x,
                top: lbl.y - lbl.height + 2.0,
                right: lbl.x + lbl.width,
                bottom: lbl.y + 2.0
            };
        }

        function anchor_boundary(a) {
            return {
                left: a.x - a.r,
                top: a.y - a.r,
                right: a.x + a.r,
                bottom: a.y + a.r
            };
        }

        function ol_area(b1, b2) {
            var x_overlap = Math.max(0, Math.min(b1.right, b2.right) - Math.max(b1.left, b2.left));
            var y_overlap = Math.max(0, Math.min(b1.bottom, b2.bottom) - Math.max(b1.top, b2.top));
            return x_overlap * y_overlap;
        }

        function energy(index) {
            // energy function, tailored for label placement
            var ener = 0,
                dx = lab[index].x - anc[index].x,
                dy = anc[index].y - lab[index].y,
                dist = Math.sqrt(dx * dx + dy * dy);

            // penalty for length of leader line
            if (dist > 0) {
                ener += dist * w_len;
            }

            // label orientation bias
            dx /= dist;
            dy /= dist;
            if (dx > 0 && dy > 0) {
                ener += 0 * w_orient;
            } else if (dx < 0 && dy > 0) {
                ener += 1 * w_orient;
            } else if (dx < 0 && dy < 0) {
                ener += 2 * w_orient;
            } else {
                ener += 3 * w_orient;
            }

            var index_boundary = label_boundary(lab[index]);

            lab.forEach(function(this_label, i) {
                if (i !== index) {
                    // penalty for intersection of leader lines
                    var overlap = intersect(anc[index], lab[index], anc[i], this_label);
                    if (overlap) {
                        ener += w_inter;
                    }

                    // penalty for label-label overlap
                    ener += (ol_area(label_boundary(this_label), index_boundary) * w_lab2);
                }

                // penalty for label-anchor overlap
                ener += (ol_area(anchor_boundary(anc[i]), index_boundary) * w_lab_anc);
            });
            return ener;
        }

        function translate(i) {
            // random translation
            lab[i].x += (Math.random() - 0.5) * max_move;
            lab[i].y += (Math.random() - 0.5) * max_move;
        }

        function rotate(i) {
            // random angle
            var angle = (Math.random() - 0.5) * max_angle;

            var s = Math.sin(angle);
            var c = Math.cos(angle);

            // translate label (relative to anchor at origin):
            lab[i].x -= anc[i].x
            lab[i].y -= anc[i].y

            // rotate label
            var lnew = { x: lab[i].x * c - lab[i].y * s, y: lab[i].x * s + lab[i].y * c };

            // translate label back
            lab[i].x = lnew.x + anc[i].x
            lab[i].y = lnew.y + anc[i].y
        }

        function mc(currT, move) {
            // Monte Carlo move

            // select a random label
            var i = Math.floor(Math.random() * lab.length);

            // save old coordinates
            var lold = { x: lab[i].x, y: lab[i].y };

            // old energy
            var old_energy = user_energy ? user_defined_energy(i, lab, anc) : energy(i);

            move(i);

            // hard wall boundaries
            if (lab[i].x > boundary.right || lab[i].x < boundary.left) {
                lab[i].x = lold.x;
            }
            if (lab[i].y > boundary.bottom || lab[i].y < boundary.top) {
                lab[i].y = lold.y;
            }

            // new energy
            var new_energy = user_energy ? user_defined_energy(i, lab, anc) : energy(i);

            if (Math.random() < Math.exp((old_energy - new_energy) / currT)) {
                acc += 1;
            } else {
                // move back to old coordinates
                lab[i].x = lold.x;
                lab[i].y = lold.y;
                rej += 1;
            }
        }


        function intersect(p1, p2, p3, p4) {
            // returns true if two lines intersect, else false
            // from http://paulbourke.net/geometry/lineline2d/
            var mua, mub;
            var denom = (p4.y - p3.y) * (p2.x - p1.x) - (p4.x - p3.x) * (p2.y - p1.y),
                numera = (p4.x - p3.x) * (p1.y - p3.y) - (p4.y - p3.y) * (p1.x - p3.x),
                numerb = (p2.x - p1.x) * (p1.y - p3.y) - (p2.y - p1.y) * (p1.x - p3.x);
            // Is the intersection along the segments
            mua = numera / denom;
            mub = numerb / denom;
            return !(mua < 0 || mua > 1 || mub < 0 || mub > 1);
        }

        function cooling_schedule(currT, initialT, nsweeps) {
            // linear cooling
            return (currT - (initialT / nsweeps));
        }

        labeler.start = function(nsweeps) {
            // main simulated annealing function
            var m = lab.length,
                currT = 1.0,
                initialT = 1.0;

            for (var i = 0; i < nsweeps; i++) {
                lab.forEach(function() {
                    var action = Math.random() < 0.5 ? translate : rotate;
                    mc(currT, action);
                });
                currT = cooling_schedule(currT, initialT, nsweeps);
            }
        };

        labeler.boundary = function(x) {
            // boundary box
            if (!arguments.length) {
                return boundary;
            }
            boundary.left = x.left;
            boundary.top = x.top;
            boundary.right = x.right;
            boundary.bottom = x.bottom;
            return labeler;
        };

        labeler.label = function(x) {
            // users insert label positions
            if (!arguments.length) {
                return lab;
            }
            lab = x;
            return labeler;
        };

        labeler.anchor = function(x) {
            // users insert anchor positions
            if (!arguments.length) {
                return anc;
            }
            anc = x;
            return labeler;
        };

        labeler.alt_energy = function(x) {
            // user defined energy
            if (!arguments.length) {
                return energy;
            }
            user_defined_energy = x;
            user_energy = true;
            return labeler;
        };

        labeler.alt_schedule = function(x) {
            // user defined cooling_schedule
            if (!arguments.length) {
                return  cooling_schedule;
            }
            user_defined_schedule = x;
            user_schedule = true;
            return labeler;
        };

        return labeler;
    };
})();

